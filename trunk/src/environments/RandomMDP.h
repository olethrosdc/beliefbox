// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_MDP_H
#define RANDOM_MDP_H

#include "DiscreteMDP.h"
#include "Environment.h"
#include "RandomNumberGenerator.h"
#include <string>
#include <vector>

/** Random MDP
    
    Creates an MDP with the given number of states and actions.
    There are three types of staes: goal states, pit states, and normal states.
    The goal and pit states always go to an absorbing terminal state.
    You can set the reward for each type of state, but they are always deterministic.
    
    The transition are random, in a way that depends on the randomness
    parameter. When randomness = 1, then all transitions are equally
    probably.  When randomness = 0, the the MDP is deterministic, with
    each action taking you to one and only one state.
 
    
 */
class RandomMDP : public DiscreteEnvironment
{
protected:
    real randomness;
    real step_value;
    real pit_value;
    real goal_value;
    RandomNumberGenerator* rng;
    bool termination;
    uint terminal_state;

public:
    RandomMDP(uint n_actions,
              uint n_states,
              real randomness,
              real step_value,
              real pit_value,
              real goal_value,
              RandomNumberGenerator* rng,
              bool termination=true);

    virtual ~RandomMDP();

    /// Generate a new MDP
    virtual DiscreteMDP* generateMDP() const;
    
    virtual DiscreteMDP* getMDP() const
    {
        return new DiscreteMDP(*mdp);
    }

    /// put the environment in its natural state
    virtual void Reset();

    /// returns true if the action succeeds
    virtual bool Act(const int& action);

    /// Remove periods from MDP
    void AperiodicityTransform(real tau) 
    {
        mdp->AperiodicityTransform(tau);
    }
    virtual real getTransitionProbability(const int& state,
                                          const int& action,
                                          const int& next_state) const
    {
        real p = mdp->getTransitionProbability(state, action, next_state);
        //printf ("# p(%d | %d, %d) = %f\n", next_state,
        //state, action, p);
        return p;
    }
    virtual real getExpectedReward(const int& state, const int& action) const
    {
        return mdp->getExpectedReward(state, action);
    }

    
protected:
    DiscreteMDP* mdp;
};

#endif
