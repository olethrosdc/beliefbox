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

#include "BayesianHierarchicalMarkovChain.h"


BayesianHierarchicalMarkovChain::BayesianHierarchicalMarkovChain(int n_states, int n_models, float prior, bool dense)
:
    BayesianMarkovChain (n_states, n_models, prior, dense),
    weight(n_models),
    log_weight(n_models),
    Lkoi(n_models, n_states),
    P_obs(n_models, n_states),
    P_next(n_states)
{
    weight[0] = 1;
    log_weight[0] = 1;
    for (int model=0; model<n_models; ++model) {
        weight[model] = prior;
        log_weight[model] = log(weight[model]);

        for (int state=0; state<n_states; state++) {
            P_obs(model, state) = 1.0 / (real) n_states;
            Lkoi(model, state) = 1.0 / (real) n_states;
        }
    }
    for (int state=0; state<n_states; state++) {
        Pr_next[state] = 1.0 / (real) n_states;
    }
    prediction=rand()%n_states;
}

BayesianHierarchicalMarkovChain::~BayesianHierarchicalMarkovChain()
{
 
#if 0
    printf("# Killing BHMC\n");
    real w = 1;
    for (int i=0; i<n_models; ++i) {
        w *= weight[i];
            printf ("%f ", w);
    }
    printf("\n");
#endif
}


/// Get the probability of the current state.
void BayesianHierarchicalMarkovChain::ObserveNextState(int state)
{
    int top_model = std::min(n_models - 1, n_observations);
    n_observations++;
    weight[0] = 1;
    log_weight[0] = 0;
        // just calculate the posterior for the current state
    Vector posterior(n_models);
    for (int model=0; model<=top_model; ++model) {
        real log_posterior = log_weight[model] + log(P_obs(model, state)) - log (Lkoi(model, state));
        
        posterior[model] = exp(log_posterior);
            //printf ("%.1f ", log_posterior);
        log_weight[model] = log_posterior;
        weight[model] = posterior[model];
    }

    for (int model=0; model<n_models; ++model) {
        //for (int model=0; model<=top_model; ++model) {
        mc[model]->ObserveNextState(state);
    }

        //printf("\n");

        /// calculate P(x_{t+1} | B_k)
        // calculate predictions for each model
    for (int model=0; model<=top_model; ++model) {
        for (int s2=0; s2<n_states; s2++) {
            P_obs(model, s2) =  mc[model]->NextStateProbability(s2);
        }
            //printf("p(%d): ", i);
        if (model == 0) {
            for (int s2=0; s2<n_states; s2++) {
                Lkoi(model, s2) = P_obs(model, s2);
            }
        } else {
            for (int s2=0; s2<n_states; s2++) {
                Lkoi(model,s2) = weight[model] * P_obs(model, s2)
                    + (1.0 - weight[model]) * Lkoi(model - 1, s2);
            }
        }
    }
    for (int s2=0; s2<n_states; s2++) {
        P_next[s2] = Lkoi(top_model, s2);
    }
    prediction = ArgMax(P_next);
}

int BayesianHierarchicalMarkovChain::predict()
{
    return prediction;
}

