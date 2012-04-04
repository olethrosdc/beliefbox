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

#ifndef DISCRETE_HIDDEN_MARKOV_MODEL_EM_H
#define DISCRETE_HIDDEN_MARKOV_MODEL_EM_H

#include "DiscreteHiddenMarkovModel.h"
#include "ExpectationMaximisation.h"
#include "RandomNumberGenerator.h"

/// This implements the standard EM algorithm with an adjustable number of
/// iterations per observation
class DiscreteHiddenMarkovModelEM
{
protected:
    int n_states;
    int n_observations;
    DiscreteHiddenMarkovModel* hmm;
    ExpectationMaximisation<DiscreteHiddenMarkovModel, int>* EM_algo;
    int n_iter;
    MultinomialDistribution B; ///< current state belief
    int T;
public:
    DiscreteHiddenMarkovModelEM(int n_states_,
                                int n_observations_,
                                real stationarity,
                                RandomNumberGenerator* rng,
                                int n_iter_ = 1) 
        : 
        n_states(n_states_),
        n_observations(n_observations_),
        B(n_states)
    {
        hmm = MakeRandomDiscreteHMM(n_states, n_observations, stationarity, rng);
        EM_algo = new ExpectationMaximisation<DiscreteHiddenMarkovModel, int>(*hmm);
        n_iter = n_iter_;
        T = 0;
    }
    ~DiscreteHiddenMarkovModelEM()
    {
        delete EM_algo;
        delete hmm;
    }
    
    
    real Observe(int x)
    {
        real ret = LOG_ZERO;
        EM_algo->Observe(x);
        T++;
        if (T > 1) {
            for (int i=0; i<1; ++i) {
                ret = EM_algo->Iterate(n_iter);
                Matrix& Ps = hmm->getBelief();
                int t = Ps.Rows();
                B = Ps.getRow(t - 1);
            }
        }
        return ret;
    }
    Vector getPrediction()
    {    
        Vector Ps(n_states);
        for (int s2=0; s2<n_states; ++s2) {
            for (int s=0; s<n_states; ++s) {
                Ps[s2] += B.Pr(s) * hmm->PrS(s,s2);
            }
        }
        int n_observations = hmm->getNObservations();
        Vector Px(n_observations);
        for (int x=0; x<n_observations; ++x) {
            Px[x] = 0;
            for (int s=0; s<n_states; ++s) {
                Px[x] += Ps(s) * hmm->PrX(s,x);
            }
        }
        return Px;
    }

    int predict()
    {
        return ArgMax(getPrediction());
    }
    void Reset()
    {
        for (int i=0; i<n_states; ++i) {
            B.Pr(i) = 1.0 / (real) n_states;
        }
        hmm->Reset();
        EM_algo->Reset();
    }

};

#endif
