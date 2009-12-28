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

BayesianPredictiveStateRepresentation::BayesianPredictiveStateRepresentation(int n_obs, int n_models, float prior, bool polya_, bool dense)
    : P_obs(n_models, n_obs), 
      Lkoi(n_models, n_obs),
      weight(n_models)
{
    beliefs.resize(n_models);
    for (int i=0; i<n_models; ++i) {
        Pr[i] = prior;
        log_prior[i] = log(Pr[i]);
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


/// Adapt the model given the next state
void BayesianPredictiveStateRepresentation::ObserveNextState(int obs)
{

//    Matrix P_obs(n_models, n_obs);
//    Matrix Lkoi(n_models, n_obs);

    int top_model = std::min(n_models - 1, n_observations);

    // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
        for (int j=0; j<n_obs; j++) {
            P_obs(model, j) =  mc[model]->NextStateProbability(j);
        }
        if (model == 0) {
            weight[model] = 1;
            for (int j=0; j<n_obs; j++) {
                Lkoi(model,j) = P_obs(model,j);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(model));
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

      
  // insert new observations
    n_observations++;
    Vector posterior(n_models);
    
    for (int model=0; model<=top_model; ++model) {
        posterior[model] = weight[model] * P_obs(model, obs) / Lkoi(model, obs);
        set_belief_param(model, log(posterior[model]) - log_prior[model]);
        mc[model]->ObserveNextObs(obs);
    }
    
}

/// Get the probability of the next state
real BayesianPredictiveStateRepresentation::NextStateProbability(int obs)
{
    Pr_next /= Pr_next.Sum();
    return Pr_next[obs];
}


/// Predict the next state
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  
///
int BayesianPredictiveStateRepresentation::predict()
{
//    Matrix P_obs(n_models, n_obs);
//    Matrix Lkoi(n_models, n_obs);
    int top_model = std::min(n_models - 1, n_observations);

        // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
        std::vector<real> p(n_obs);

        for (int obs=0; obs<n_obs; obs++) {
            P_obs(model, obs) =  mc[model]->NextObsProbability(obs);
        }
            //printf("p(%d): ", i);
        if (model == 0) {
            weight[model] = 1;
            for (int obs=0; obs<n_obs; obs++) {
                Lkoi(model, obs) = P_obs(model, obs);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(model));
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


