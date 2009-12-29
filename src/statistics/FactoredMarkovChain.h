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

/** A factored Markov chain.

    This is not really a statistical distribution as it is composed of
    two variables, one of which is externally controllable.

    The variable \f$x = (y, z)\f$ is a pair with the controllable
    action \f$y\f$ arising from process \f$\pi\f$ and the random
    observation \f$z\f$ being a predicted variable arising from
    process \f$\mu\f$.
    
    Thus, the probability distribution over next observations
    $x_{t+1}$ given the history $x^t$ is can be written as:

    \f[
    \Pr(x_{t+1} | x^t, \mu, \pi) = 
    \Pr(y_{t+1} | z_{t+1}, x^t, \pi)
    \Pr(z_{t+1} | x^t, \mu)
    \f]

    Thus, the next observation depends on the complete history of
    observations and actions.
 */
class FactoredMarkovChain 
{
protected:
    int n_actions; ///< number of actions
    int n_obs; ///< number of observations
    int mem_size; ///< order of the chain
	SparseTransitions transitions; ///< history-wide transition table
    std::vector<int> act_history;
    std::vector<int> obs_history;
    std::vector<int> history;
    real threshold;
public:
    FactoredMarkovChain(int n_actions_, int n_obs_, int mem_size);
    virtual ~FactoredMarkovChain();

    /* probabilities */
    real getTransition(FactoredContext ctx, int prd);
    real getProbability(FactoredContext ctx, int prd);
    void getProbabilities(FactoredContext ctx, std::vector<real>& p);
    void getNextStateProbabilities(int act, std::vector<real>& p);

    /* Training and generation */
    int PushObservation (int act, int prd);
    real Observe (int act, int prd);
    real NextStateProbability (int act, intstate);
    void Reset();
    int GenerateStatic();
    int GenerateStatic(int act);
    
    /* Debug functions */
    virtual int ShowTransitions ();
}; 
/*@}*/
#endif
