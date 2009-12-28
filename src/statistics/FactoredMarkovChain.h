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
#ifndef FACTORED_MARKOVCHAIN_H
#define FACTORED_MARKOVCHAIN_H

#include <map>
#include <vector>
#include "real.h"
#include "SparseTransitions.h"
/**
   \ingroup StatisticsGroup
 */
/*@{*/

typedef long FactoredContext;

/** A factored Markov chain

    The observable \f$x = (y, z)\f$ is a pair with the controllable
    action \f$y\f$ element being a controllable variable arising from
    process \f$\pi\f$ and the random \f$z\f$ element being a predicted
    variable arising from process \f$\mu\f$.
    
    Thus, the probability distribution over next observations
    $x_{t+1}$ given the history $x^t$ is can be written as:

    \f[
    \Pr(x_{t+1} | x^t, \mu, \pi) = 
    \Pr(y_{t+1} | x^t, \pi)
    \Pr(z_{t+1} | y_{t+1}, x^t, \mu)
    \f]
 */
class FactoredMarkovChain 
{
protected:
    int n_transitions; ///< total number of transitions
	SparseTransitions transitions; ///< history-wide transition table

    real threshold;
public:
    FactoredMarkovChain (int n_states, int mem_size);
    virtual ~FactoredMarkovChain ();

    /* probabilities */
    real getTransition (FactoredContext ctx, int prd);
    real getProbability (FactoredContext ctx, int prd);
    void getProbabilities(FactoredContext ctx, std::vector<real>& p);
    void getNextStateProbabilities(std::vector<real>& p);

    /* Training and generation */
    int PushState (int prd);
    real ObserveNextState (int state);
    real NextStateProbability (int state);
    void Reset();
    int GenerateStatic();
    int generate();
    
    /* Debug functions */
    virtual int ShowTransitions ();
}; 
/*@}*/
#endif
