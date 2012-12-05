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
#include "Vector.h"
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
	bool endsim; ///< Absorbing state
    uint n_states; ///< The state dimension
    uint n_actions; ///< The action dimension
    S state_lower_bound; ///< lower bound on the states
    S state_upper_bound; ///< upper bound on the states
public:
    Environment() : n_states(1), n_actions(1)
    {
        state_lower_bound = 0;
        state_upper_bound = 0;
        reward = 0.0;
		endsim = false;
    }

    Environment(int n_states_, int n_actions_)
  : n_states(n_states_), n_actions(n_actions_)
    {
        state_lower_bound = 0;
        state_upper_bound = n_states;
        reward = 0.0;
		endsim = false;
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
    virtual bool Act(const A& action) = 0;
    
    /// Return a full MDP model of the environment.  This may not be
    /// possible for some environments.  The MDP is required to be
    /// freed by the user!
    virtual MDP<S, A>* getMDP() const
    {
        return NULL;
    }
    virtual const char* Name() const
    {
        return "Undefined environment name";
    }
	// --- The following functions are not supposed to be overwritten.. -- //
    /// returns a (reference to) the current state
    const S& getState() const
    {
        return state;
    }
	///  sets the current state
	void setState(const S& s_next)
	{
		state = s_next;
	}
    /// returns the current reward
    real getReward() const
    {
        return reward;
    }
	/// indicates if the current state is absorbing or not
	bool getEndsim() const
	{
		return endsim;
	}
    /// returns the number of state dimensions
    uint getNStates() const
    {
        return n_states;
    }
    uint getNActions() const
    {
        return n_actions;
    }
    /// Set the overall randomness of the environment
    virtual void setRandomness(real randomness)
    {
        
    }
    const S& StateUpperBound() const
    {
        return state_upper_bound;
    }
    const S& StateLowerBound() const
    {
        return state_lower_bound;
    }
	virtual real getTransitionProbability(const S& state, const A& action, const S& next_state) const = 0;

    virtual real getExpectedReward(const S& state, const A& action) const = 0;

};

/// Default type for discrete environments
typedef Environment<int, int> DiscreteEnvironment;

/// Default type for continuous state environments
typedef Environment<Vector, int> ContinuousStateEnvironment;

/// Default type for continuous environments
typedef Environment<Vector, Vector> ContinuousEnvironment;

/*@}*/

#endif
