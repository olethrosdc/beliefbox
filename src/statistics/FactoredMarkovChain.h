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
#include <cassert>
#include "real.h"
#include "Ring.h"
#include "SparseTransitions.h"
#include "FactoredPredictor.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/



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
class FactoredMarkovChain  : public FactoredPredictor
{
public:
    typedef long long Context;
protected:
    int T; ///< total time passed
    int n_actions; ///< number of actions
    int n_obs; ///< number of observations
    int n_states; ///< number of action*observations
    int mem_size; ///< order of the chain
    Context n_contexts; ///< number of possible total contexts
    SparseTransitions transitions; ///< history-wide transition table

    /// The current context. It is updated whenever a new action
    /// \f$a_t\f$ is observed and is a summary of the history of the
    /// past \f$(x_{t-i}, a_{t-i})_{i=0}^{\tau-1}\f$ observation pairs.
    Context current_context;
    Ring<int> act_history;
    Ring<int> obs_history;
    Ring<int> history;
    real threshold;

    /// Calculate the local state for a given action a, observation x
    int CalculateState(int a, int x)
    {
        assert(a>=0 && a<n_actions);
        assert(x>=0 && x<n_obs);
        return a*n_obs + x;
    }

    /// a state is pushed only at the end
    void PushState(int state)
    {
        assert(state>=0 && state<n_states);
        history.push_back(state);
        T++;
        assert(T==history.size() && T==act_history.size() && T==obs_history.size());
    }

    /// You push an action first, and then the observation
    void PushAction(int act)
    {
        act_history.push_back(act);
        PushState(CalculateState(act, obs_history.back()));
        assert(act_history.size()==obs_history.size());
    }

    void PushObservation(int obs)
    {
        assert(act_history.size()==obs_history.size());
        obs_history.push_back(obs);
        //PushState(CalculateState(act, obs_history.back()));
    }
public:
    FactoredMarkovChain(int n_actions_, int n_obs_, int mem_size);
    virtual ~FactoredMarkovChain();

    /* house keeping */
    Context CalculateContext();
    Context getContext(int act);

    /* probabilities */
    real getTransition(Context context, int prd);
    real getProbability(Context context, int prd);
    void getProbabilities(Context context, std::vector<real>& p);
    void getNextStateProbabilities(int act, std::vector<real>& p);
    
    /* Training and generation */
  virtual real Observe (int prd);
  virtual real Observe (int act, int prd);
  virtual real ObservationProbability (int act, int x);
  virtual real ObservationProbability (int x);
  virtual void Reset();
    int GenerateStatic();
    int GenerateStatic(int act);
    
}; 
/*@}*/
#endif
