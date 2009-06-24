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

#include "BayesianMarkovChain.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "Random.h"

BayesianMarkovChain::BayesianMarkovChain(int n_states, int n_models, float prior, bool dense)
{
    this->n_models = n_models;
    this->n_states = n_states;

    mc.resize(n_models);
    log_prior.resize(n_models);

    Pr.Resize(n_models);
    Pr_next.Resize(n_states);

    float sum = 0.0;
    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = 1.0 / (float) n_states;
    }

    for (int i=0; i<n_models; ++i) {
        Pr[i] = 1.0 / (1.0 + (float) i);
        log_prior[i] = log(Pr[i]);
        sum += Pr[i];
        if (dense) {
            mc[i] = new DenseMarkovChain(n_states, i);
        } else {
            mc[i] = new SparseMarkovChain(n_states, i);
        }
        mc[i]->setThreshold(prior);
    }
    for (int i=0; i<n_models; ++i) {
        Pr[i] /= sum;
    }
}

BayesianMarkovChain::~BayesianMarkovChain()
{
    for (int i=0; i<n_models; ++i) {
        delete mc[i];
    }
}

void BayesianMarkovChain::Reset()
{
    for (int i=0; i<n_models; ++i) {
        mc[i]->Reset();
    }
}

/// Get the probability of the current state.
void BayesianMarkovChain::ObserveNextState(int state)
{
    float sum = 0.0;
    for (int i=0; i<n_models; ++i) {
        Pr[i] *= mc[i]->NextStateProbability(state);
        sum += Pr[i];
    }

    for (int i=0; i<n_models; ++i) {
        mc[i]->ObserveNextState(state);
        Pr[i] /= sum;
    }
}

/// Get the probability of the next state
float BayesianMarkovChain::NextStateProbability(int state)
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
/// samples (does it ?)
///
/// Side-effects: Changes the current state.
int BayesianMarkovChain::generate()
{
    int i = generate_static();
    for (int j=0; j<n_models; ++j) {
        mc[j]->PushState(i);
    }
    return i;
}


/// Generate the next state.
///
/// We are flattening the hierarchical distribution to a simple
/// multinomial.  This allows us to more accurately generate random
/// samples (does it ?)
///
int BayesianMarkovChain::generate_static()
{
    float sum = 0.0;
    for (int i=0; i<n_states; ++i) {
        Pr_next[i] = 0.0;
        for (int j=0; j<n_models; ++j) {
            Pr_next[i] += Pr[j] * mc[j]->NextStateProbability(i);
        }
    }

   sum = 0.0;
    float X = urandom();
    for (int i=0; i<n_states; ++i) {
        sum += Pr_next[i];
        if (X<sum && Pr_next[i] > 0.0f) {
            return i;
        }
    }

    //Swarning ("Multinomial generation failed.");
    int i = rand()%n_states;
    return i;
}

void BayesianMarkovChain::ShowTransitions()
{
    for (int j=0; j<n_models; ++j) {
        mc[j]->ShowTransitions();
    }
}
