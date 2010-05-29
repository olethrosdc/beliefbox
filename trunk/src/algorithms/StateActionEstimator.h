/* -*- Mode: C++; -*- */
/* VER: $Id: StateActionEstimator.h,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef STATE_ACTION_ESTIMATOR_H
#define STATE_ACTION_ESTIMATOR_H

#include "real.h"

class StateActionEstimator
{
 public:
    int n_states;  ///< number of states
    int n_actions; ///< number of actions

    StateActionEstimator(int n_states, int n_actions)
    {
        this->n_states = n_states;
        this->n_actions = n_actions;
    }
    /// Destructor
    virtual ~StateActionEstimator()
    {
    }
    /// Observe a state-action-reward-state tuplet
    /// s_p : previous state
    /// a : action
    /// r : reward
    /// s : state
    virtual void Observe(int s_p, int a_p, real r, int s, int a) = 0;

    /// Get transition probability.
    /// When we have a distribution over transition probabilities then
    /// it returns the mean of the distribution.
    virtual real TransProbability(int s_p, int a_p, int s) = 0;

    /// Sample from the transition distribution.
    /// Again, here it samples from the mean when we have a prior.
    virtual int SampleTransition(int s, int a) = 0;

    /// If (s_p=i, a=j, s=k), then get the density of the ijk 
    /// transition being x under our prior \f$\xi\f$:
    /// \f[
    /// p( P(ijk)=x | \xi)
    /// \f]
    virtual real GetTransitionPrior(int s_p, int a_p, int s, real x) = 0;
    
    /// Sample a possible value for the transition probability
    /// according to our current prior.
    virtual real SampleTransitionPrior(int s_p, int a_p, int s) = 0;

    /// Reset the estimator
    virtual void Reset() = 0;
};

#endif
