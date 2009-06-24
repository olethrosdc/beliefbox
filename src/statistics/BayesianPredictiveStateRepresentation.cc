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

BayesianPredictiveStateRepresentation::BayesianPredictiveStateRepresentation(int n_states, int n_models, float prior, bool dense)
    : BayesianMarkovChain(n_states, n_models, prior, dense)
{
    beliefs.resize(n_models);
}

BayesianPredictiveStateRepresentation::~BayesianPredictiveStateRepresentation()
{
}

void BayesianPredictiveStateRepresentation::Reset()
{
    for (int i=0; i<n_models; ++i) {
        mc[i]->Reset();
    }
}


/// Get the probability of the current state.
void BayesianPredictiveStateRepresentation::ObserveNextState(int state)
{
    std::vector<real> weight(n_models);
    Matrix Lkoi(n_models, n_states);

    for (int i=0; i<n_models; ++i) {
        std::vector<real> p(n_states);
        mc[i]->getNextStateProbabilities(p);            
        if (i == 0) {
            weight[i] = 1;
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = p[j];
            }
        } else {
            int state_i = mc[i]->getCurrentState();
            weight[i] = exp(log_prior[i] + get_belief_param(i, state_i));
            for (int j=0; j<n_states; j++) {
                Lkoi(i,j) = weight[i] * p[j] + (1.0 - weight[i])*Lkoi(i,j-1);
            }
        }
    }

    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = 0.0;
        for (int j=0; j<n_models; ++j) {
            Pr_next[i] += Pr[j] * mc[j]->NextStateProbability(i);
        }
    }
}

/// Get the probability of the next state
float BayesianPredictiveStateRepresentation::NextStateProbability(int state)
{
    float sum = 0.0;
    for (int i=0; i<n_models; ++i) {
        sum += Pr[i] * mc[i]->NextStateProbability(state);
    }
    return sum;
}

/// Generate the next state.
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  This allows us to more accurately generate random
/// samples (!) does it ?
///
int BayesianPredictiveStateRepresentation::generate_static()
{
    float sum = 0.0;
    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = 0.0;
        for (int j=0; j<n_models; ++j) {
            Pr_next[i] += Pr[j] * mc[j]->NextStateProbability(i);
        }
        //sum  += Pr_next[i];
    }

   sum = 0.0;
    float X = urandom();
    for (int i=0; i<n_states; ++i) {
        sum += Pr_next[i];
        if (X<sum && Pr_next[i] > 0.0f) {
            for (int j=0; j<n_models; ++j) {
                mc[j]->PushState(i);
            }
            return i;
        }
    }

    //Swarning ("Multinomial generation failed.");
    int i = rand()%n_states;
    for (int j=0; j<n_models; ++j) {
        mc[j]->PushState(i);
    }
    return i;
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
    int i = generate_static();
    for (int j=0; j<n_models; ++j) {
        mc[j]->PushState(i);
    }
    return i;
}

