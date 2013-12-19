// -*- Mode: c++ -*-
// copyright (c) 2007-2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TRANSITION_DISTRIBUTION_H
#define TRANSITION_DISTRIBUTION_H

#include "real.h"
#include <map>

template <typename StateType, typename ActionType>
class TransitionDistribution
{
public:
    virtual ~TransitionDistribution() {}
    virtual StateType generate(StateType state, ActionType action) = 0;
    virtual real pdf(StateType state, ActionType action, StateType next_state) = 0;
};

template <typename StateType, typename ActionType>
struct Transition
{
	StateType state;
	ActionType action;
	StateType next_state;
	Transition(StateType s, 
			   ActionType a,
			   StateType s2)
		: state(s), action(a), next_state(s2)
	{}
};


 
/** Discrete transition distribution.

	In this model, there is no remaining probability mass. So, the probability of going to any particular state is zero.
	
 */
template<>
class TransitionDistribution<int, int>
{
public:
	int n_states; ///< the maximum number of state
	int n_actions; ///< the maximum number of actions
	std::map<Transition<int, int>, real> P;
	TransitionDistribution(int n_states_, int n_actions_)
		: n_states(n_states_),
		  n_actions(n_actions_)
	{
	}
	virtual ~TransitionDistribution() {}
	/// Set a state transition
	virtual void SetTransition(int state, int action, int next_state, real probability);
	/// Generate a next state
	virtual int generate(int state, int action);
	/// Get the probability of the next state
	virtual real pdf(int state, int action, int next_state);

};

#endif
