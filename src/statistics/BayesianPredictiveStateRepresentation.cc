/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BayesianPredictiveStateRepresentation.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "Random.h"
#include "Matrix.h"
#include "Distribution.h"

BayesianPredictiveStateRepresentation::BayesianPredictiveStateRepresentation(int n_obs_, int n_actions_, int n_models_, float prior)
    : total_observations(0),
      n_obs(n_obs_),
      n_actions(n_actions_),
      n_models(n_models_),
      P_obs(n_models, n_obs), 
      Lkoi(n_models, n_obs),
      weight(n_models),
      mc(n_models),
      log_prior(n_models),
      Pr(n_models),
      Pr_next(n_obs)
{
    beliefs.resize(n_models);
    for (int i=0; i<n_models; ++i) {
        Pr[i] = prior;
        log_prior[i] = log(Pr[i]);
        mc[i] = new FactoredMarkovChain(n_actions, n_obs, i);
    }
}

BayesianPredictiveStateRepresentation::~BayesianPredictiveStateRepresentation()
{
    //printf("Killing BPSR\n");
}

void BayesianPredictiveStateRepresentation::Reset()
{
    for (int i=0; i<n_models; ++i) {
        mc[i]->Reset();
    }
}


/** Seed the model with an initial observation.
 
    @param observation the observation
*/
void BayesianPredictiveStateRepresentation::Observe(int observation)
{
    for (int model=0; model<n_models; ++model) {
        mc[model]->Observe(observation);
    }
    
}

/** Adapt the model given a specific action and next observation.
 
    @param action the action
    @param observation the observation
*/
void BayesianPredictiveStateRepresentation::Observe(int action, int observation)
{

    //    Matrix P_obs(n_models, n_obs);
    //    Matrix Lkoi(n_models, n_obs);

    int top_model = std::min(n_models - 1, total_observations);
    
    // calculate predictions for each model for the given action 
    for (int model=0; model<=top_model; ++model) {
        for (int j=0; j<n_obs; j++) {
            P_obs(model, j) =  mc[model]->ObservationProbability(action, j);
        }
        if (model == 0) {
            weight[model] = 1;
            for (int j=0; j<n_obs; j++) {
                Lkoi(model,j) = P_obs(model,j);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(action, model));
            for (int j=0; j<n_obs; j++) {
                Lkoi(model,j) = weight[model] * P_obs(model, j) + (1.0 - weight[model])*Lkoi(model-1, j); 
            }
        }
    }
    real p_w = 1.0;
    for (int model=top_model; model>=0; model--) {
        Pr[model] = p_w * weight[model];
        p_w *= (1.0 - weight[model]);

    }

    real sum_pr_s = 0.0;
    for (int s=0; s<n_obs; ++s) {
        real Pr_s = 0;
        for (int model=0; model<=top_model; ++model) {
            Pr_s += Pr[model]*P_obs(model, s);
        }
        Pr_next[s] = Pr_s;
        sum_pr_s += Pr_s;
    }

      
    // insert new observations to all models
    total_observations++;
    Vector posterior(n_models);
    
    for (int model=0; model<=top_model; ++model) {
        posterior[model] = weight[model] * P_obs(model, observation) / Lkoi(model, observation);
        set_belief_param(action, model, log(posterior[model]) - log_prior[model]);
        mc[model]->Observe(action, observation); ///< NOTE: Maybe this should be in a different loop?
    }
    
}

/** Get the probability of the next state
    
    P_t(x_{t+1} = x | a_t = a)

*/
real BayesianPredictiveStateRepresentation::ObservationProbability(int action, int observation)
{
    int top_model = std::min(n_models - 1, total_observations);
    printf("models:%d - top: %d\n", n_models, top_model);
    // calculate predictions for each model for the given action 
    for (int model=0; model<=top_model; ++model) {
        for (int j=0; j<n_obs; j++) {
            P_obs(model, j) =  mc[model]->ObservationProbability(action, j);
            //printf ("%d %d %f\n", model, j, P_obs(model, j));
        }
        if (model == 0) {
            weight[model] = 1;
            for (int j=0; j<n_obs; j++) {
                Lkoi(model,j) = P_obs(model,j);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(action, model));
            for (int j=0; j<n_obs; j++) {
                Lkoi(model,j) = weight[model] * P_obs(model, j) + (1.0 - weight[model])*Lkoi(model-1, j); 
            }
        }
    }
    real p_w = 1.0;
    for (int model=top_model; model>=0; model--) {
        Pr[model] = p_w * weight[model];
        p_w *= (1.0 - weight[model]);

    }

    real sum_pr_s = 0.0;
    for (int s=0; s<n_obs; ++s) {
        real Pr_s = 0;
        for (int model=0; model<=top_model; ++model) {
		  //printf ("%d %d %f\n", model, s, P_obs(model, s));
            Pr_s += Pr[model]*P_obs(model, s);
        }
        Pr_next[s] = Pr_s;
        sum_pr_s += Pr_s;
    }

    return Pr_next[observation];
}


/// Predict the next state
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  
///
int BayesianPredictiveStateRepresentation::predict(int a)
{
    //    Matrix P_obs(n_models, n_obs);
    //    Matrix Lkoi(n_models, n_obs);
    int top_model = std::min(n_models - 1, total_observations);

    // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
        std::vector<real> p(n_obs);

        for (int obs=0; obs<n_obs; obs++) {
            P_obs(model, obs) =  mc[model]->ObservationProbability(a, obs);
        }
        //printf("p(%d): ", i);
        if (model == 0) {
            weight[model] = 1;
            for (int obs=0; obs<n_obs; obs++) {
                Lkoi(model, obs) = P_obs(model, obs);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(a, model));
            for (int obs=0; obs<n_obs; obs++) {
                Lkoi(model,obs) = weight[model] * P_obs(model, obs)
                    + (1.0 - weight[model]) * Lkoi(model - 1, obs);
            }
        }
    }

    real p_w = 1.0;
    for (int model=top_model; model>=0; model--) {
        Pr[model] = p_w * weight[model];
        p_w *= (1.0 - weight[model]);
    }
    
    Pr /= Pr.Sum(0, top_model);
    for (int obs=0; obs<n_obs; ++obs) {
        Pr_next[obs] = 0.0;
        for (int model=0; model<=top_model; model++) {
            Pr_next[obs] += Pr[model]*P_obs(model, obs);
        }
    }

    return ArgMax(&Pr_next);
    //return DiscreteDistribution::generate(Pr_next);
}


