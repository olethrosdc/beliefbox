/* -*- Mode: c++;  -*- */
/*VER: $Id: MarkovChain.h,v 1.7 2006/08/17 05:35:00 olethros Exp $*/
// copyright (c) 2004 by Christos Dimitrakakis <dimitrak@idiap.ch>
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

BayesianPredictiveStateRepresentation::BayesianPredictiveStateRepresentation(int n_states, int n_models, float prior, bool polya_, bool dense)
    : BayesianMarkovChain(n_states, n_models, prior, dense),
      polya(polya_),
      P_obs(n_models, n_states), 
      Lkoi(n_models, n_states),
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
void BayesianPredictiveStateRepresentation::ObserveNextState(int state)
{

//    Matrix P_obs(n_models, n_states);
//    Matrix Lkoi(n_models, n_states);

    int top_model = std::min(n_models - 1, n_observations);

#if 1
        // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
            //std::vector<real> p(n_states);

        for (int j=0; j<n_states; j++) {
            P_obs(model, j) =  mc[model]->NextStateProbability(j);
        }
            //printf("p(%d): ", i);
        if (model == 0) {
            weight[model] = 1;
            for (int j=0; j<n_states; j++) {
                Lkoi(model,j) = P_obs(model,j);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(model));
            for (int j=0; j<n_states; j++) {
                Lkoi(model,j) = weight[model] * P_obs(model, j) + (1.0 - weight[model])*Lkoi(model-1, j); 
            }
        }
    }
#else
    int j = state;
    for (int model=0; model<=top_model; ++model) {
        P_obs(model, j) =  mc[model]->NextStateProbability(j);
            //printf("p(%d): ", i);
        if (model == 0) {
            weight[model] = 1;
            Lkoi(model,j) = P_obs(model,j);
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(model));
            Lkoi(model,j) = weight[model] * P_obs(model, j) + (1.0 - weight[model])*Lkoi(model-1, j); 
        }
    }
#endif
    
    real p_w = 1.0;
    for (int model=top_model; model>=0; model--) {
        Pr[model] = p_w * weight[model];
        p_w *= (1.0 - weight[model]);

    }

    real sum_pr_s = 0.0;
    for (int s=0; s<n_states; ++s) {
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
        if (polya) {
            real par = exp(get_belief_param(model));// + log_prior[model]);
            set_belief_param(model, log(1.0 + par));// - log_prior[model]);
        } else {
            posterior[model] = weight[model] * P_obs(model, state) / Lkoi(model, state);
            set_belief_param(model, log(posterior[model]) - log_prior[model]);
        }
        mc[model]->ObserveNextState(state);
    }
    
}

/// Get the probability of the next state
real BayesianPredictiveStateRepresentation::NextStateProbability(int state)
{
    Pr_next /= Pr_next.Sum();
    return Pr_next[state];
}


/// Predict the next state
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  
///
int BayesianPredictiveStateRepresentation::predict()
{
//    Matrix P_obs(n_models, n_states);
//    Matrix Lkoi(n_models, n_states);
    int top_model = std::min(n_models - 1, n_observations);

        // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
        std::vector<real> p(n_states);

        for (int state=0; state<n_states; state++) {
            P_obs(model, state) =  mc[model]->NextStateProbability(state);
        }
            //printf("p(%d): ", i);
        if (model == 0) {
            weight[model] = 1;
            for (int state=0; state<n_states; state++) {
                Lkoi(model, state) = P_obs(model, state);
            }
        } else {
            weight[model] = exp(log_prior[model] + get_belief_param(model));
            for (int state=0; state<n_states; state++) {
                Lkoi(model,state) = weight[model] * P_obs(model, state)
                    + (1.0 - weight[model]) * Lkoi(model - 1, state);
            }
        }
    }

    real p_w = 1.0;
    for (int model=top_model; model>=0; model--) {
        Pr[model] = p_w * weight[model];
        p_w *= (1.0 - weight[model]);
    }
#if 0
    for (int i=0; i<=top_model; i++) {
        printf ("%f ", Pr[i]);
    }
    printf("#BPSR \n");
#endif
    
    Pr /= Pr.Sum(0, top_model);
    for (int state=0; state<n_states; ++state) {
        Pr_next[state] = 0.0;
        for (int model=0; model<=top_model; model++) {
            Pr_next[state] += Pr[model]*P_obs(model, state);
        }
    }

#if 0
    for (int model=0; model<=top_model; model++) {
        real s = P_obs.RowSum(model);
        if (fabs(s - 1.0) > 0.001) {
            if (polya) {
                fprintf (stderr, "polya ");
            }
            fprintf(stderr, "sum[%d]: %f\n", model, s);
        }
    }
    real sum_pr_s = Pr_next.Sum();
    if (fabs(sum_pr_s - 1.0) > 0.001) {
        if (polya) {
            fprintf (stderr, "polya ");
        }
        fprintf(stderr, "sum2: %f\n", sum_pr_s);
        fprintf(stderr, "model sum: %f\n", Pr.Sum(0, top_model));
        //exit(-1);
    }
#endif
 
    
#if 0
    printf ("# BPSR ");
    for (int i=0; i<n_states; ++i) {
        printf ("(%f %f", Pr_next[i], Lkoi(top_model,i));
    }
    printf("\n");
#endif
    return ArgMax(&Pr_next);
    //return DiscreteDistribution::generate(Pr_next);
}

/// Generate the next state.
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  This allows us to more accurately generate random
/// samples (!) does it ?
///
/// Side-effects: Changes the current state.
int BayesianPredictiveStateRepresentation::generate()
{
    int i = predict();
    for (int j=0; j<n_models; ++j) {
        mc[j]->PushState(i);
    }
    return i;
}

