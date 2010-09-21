// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "MDP.h"
#include <cstdlib>
/**
   \defgroup EnvironmentGroup Environments
 */
/**
   \ingroup EnvironmentGroup
 */
/*@{*/
/// Template for environments
template <typename S, typename A>
class Environment
{
protected:
    S state; ///< The current state
    real reward; ///< The current reward
    uint n_states; ///< The state dimension
    uint n_actions; ///< The action dimension
public:
    Environment() : n_states(1), n_actions(1)
    {
    }

    Environment(int n_states_, int n_actions_)
  : n_states(n_states_), n_actions(n_actions_)
    {
    }

    virtual ~Environment() 
    {
    }

    /// put the environment in its "natural: state
    virtual void Reset() = 0;

    /// returns true if the action succeeds, false if it does not.
    ///
    /// The usual of false is that the environment is in a terminal
    /// absorbing state.
    virtual bool Act(A action) = 0;
    
    /// Return a full MDP model of the environment. 
    /// This may not be possible for some environments
    virtual MDP<S, A>* getMDP() const
    {
        return NULL;
    }
    /// returns a (reference to) the current state
    S& getState()
    {
        return state;
    }
	///  sets the current state
	void setState(S s_next)
	{
		state = s_next;
	}
    /// returns the current reward
    real getReward()
    {
        return reward;
    }
    /// returns the number of state dimensions
    uint getNStates()
    {
        return n_states;
    }
    uint getNActions()
    {
        return n_actions;
    }
    virtual const char* Name()
    {
        return "Undefined environment name";
    }

};

/// Default type for discrete environments
typedef Environment<int, int> DiscreteEnvironment;

/*@}*/

#endif
