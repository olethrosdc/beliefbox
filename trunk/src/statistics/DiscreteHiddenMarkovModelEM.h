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
    DiscreteHiddenMarkovModel* hmm;
    ExpectationMaximisation<DiscreteHiddenMarkovModelEM, int>* EM_algo;
public:
    DiscreteHiddenMarkovModelEM(int n_states,
                                int n_observations,
                                real stationarity,
                                RandomNumberGenerator rng) 
    {
        hmm = MakeRandomDiscreteHMM(n_states, n_observations, stationarity, rng);
        EM_algo = new ExpectationMaximisation<DiscreteHiddenMarkovModel, int>(*hmm);
    }
    ~DiscreteHiddenMarkovModelEM()
    {
        delete EM_algo;
        delete hmm;
    }
};

#endif
