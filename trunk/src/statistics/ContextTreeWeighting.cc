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

#include "ContextTreeWeighting.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "Random.h"
#include "Matrix.h"
#include "Distribution.h"

ContextTreeWeighting::ContextTreeWeighting(int n_states, int n_models, float prior, bool dense)
    : BayesianMarkovChain(n_states, n_models, prior, dense)
{
    n_observations = 0;
}

ContextTreeWeighting::~ContextTreeWeighting()
{
    //printf("Killing BPSR\n");
}

void ContextTreeWeighting::Reset()
{
    for (int i=0; i<n_models; ++i) {
        mc[i]->Reset();
    }
}


/// Get the probability of the current state.
void ContextTreeWeighting::ObserveNextState(int state)
{
    std::vector<real> weight(n_models);
    Matrix Lkoi(n_models, n_states); // p obs given all models to k
    Matrix P_obs(n_models, n_states); // probability of obs for model k
    int top_model = std::min(n_models - 1, n_observations);

        // calculate predictions for each model
    for (int i=0; i<=top_model; ++i) {
        std::vector<real> p(n_states);

        for (int j=0; j<n_states; j++) {
            P_obs(i,j) =  mc[i]->NextStateProbability(j);
        }
            //printf("p(%d): ", i);
        if (i == 0) {
            weight[i] = exp(log_prior[i]);
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = P_obs(i,j);
            }
        } else {
            weight[i] = exp(log_prior[i]);
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = weight[i] * P_obs(i,j) + (1.0 - weight[i])*Lkoi(i-1,j); // NOTE: was i, j-!!!
            }
        }
    }
    
    real p_w = 1.0;
    for (int i=top_model; i>=0; i--) {
        Pr[i] = p_w*weight[i];
        p_w *= (1.0 - weight[i]);

    }

    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = P_obs(top_model, i);
        real Pr_i = 0;
        for (int j=0; j<=top_model; ++j) {
            Pr_i += Pr[j]*P_obs(j, i);
        }
        Pr_next[i] = Pr_i;
    }

        // insert new observations
    n_observations++;
    Vector posterior(n_models);
    
    for (int i=0; i<=top_model; ++i) {
        posterior[i] = weight[i] * P_obs(i, state) / Lkoi(i, state);
        mc[i]->ObserveNextState(state);
    }
    
}

/// Get the probability of the next state
float ContextTreeWeighting::NextStateProbability(int state)
{
    return Pr_next[state];
}

/// Generate the next state.,
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  This allows us to more accurately generate random
/// samples (!) does it ?
///
int ContextTreeWeighting::predict()
{
    std::vector<real> weight(n_models);
    Matrix Lkoi(n_models, n_states); // p obs given all models to k
    Matrix P_obs(n_models, n_states); // probability of obs for model k
    int top_model = std::min(n_models - 1, n_observations);

        // calculate predictions for each model
    for (int i=0; i<=top_model; ++i) {
        std::vector<real> p(n_states);

        for (int j=0; j<n_states; j++) {
            P_obs(i,j) =  mc[i]->NextStateProbability(j);
        }
            //printf("p(%d): ", i);
        if (i == 0) {
            weight[i] = 1;
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = P_obs(i,j);
            }
        } else {
            weight[i] = exp(log_prior[i]);
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = weight[i] * P_obs(i,j) + (1.0 - weight[i])*Lkoi(i-1,j); //!!!!!!!!!! NOTE: was i, j-1
            }
        }
    }

    real p_w = 1.0;
    for (int i=top_model; i>=0; i--) {
        Pr[i] = p_w*weight[i];
        p_w *= (1.0 - weight[i]);
    }
#if 0
    for (int i=0; i<n_models; i++) {
        printf ("%f ", Pr[i]);
    }
    printf("#BPSR \n");
#endif
    
    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = 0.0;
        for (int j=0; j<=top_model; j++) {
            Pr_next[i] += Pr[j]*P_obs(j, i);
        }
    }

#if 0
    for (int i=0; i<n_states; ++i) {
        printf ("%f ", Pr_next[i]);
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
int ContextTreeWeighting::generate()
{
    int i = predict();
    for (int j=0; j<n_models; ++j) {
        mc[j]->PushState(i);
    }
    return i;
}

