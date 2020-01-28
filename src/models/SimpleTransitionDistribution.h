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

#ifndef SIMPLE_TRANSITION_DISTRIBUTION_H
#define SIMPLE_TRANSITION_DISTRIBUTION_H

#include "real.h"
#include "DiscreteStateSet.h"
#include "StateAction.h"
#include <cstdio>



/** Simple Discrete transition distribution.

	In this model, we employ an unorder map of actual transitions, as well as a map of next states.
 */

class SimpleTransitionDistribution
{
public:
	int n_states; ///< the maximum number of state
	int n_actions; ///< the maximum number of actions
	DiscreteStateSet empty_set; ///< included for convenience
	/// The implementation of the discrete transition distribution
    real *arr; ///< gives the actual probabilities
	std::unordered_map<DiscreteStateAction, DiscreteStateSet> next_states; ///< next states for quick access
	SimpleTransitionDistribution(int n_states_, int n_actions_)
		: n_states(n_states_),
		  n_actions(n_actions_)
	{

    	arr = (real *)calloc(n_states*n_actions * n_states , sizeof(real));
/*
		for (int s=0; s<n_states; s++)
			for (int a=0; a<n_actions; a++)
				for (int s2=0; s2<n_states; s2++) arr[(s+a*n_states)*n_states + s2]=0.0f;
*/
	}

	virtual ~SimpleTransitionDistribution();

	int GetNStates() const
	{
		return n_states;
	}

	int GetNActions() const
	{
		return n_actions;
	}

	/// Set a state transition
	virtual void SetTransition(int state, int action, int next_state, real probability);

	/// Get a state transition
	virtual real GetTransition(int state, int action, int next_state) const;

	/// Generate a next state
	virtual int generate(int state, int action) const;
	/// Get the probability of the next state
	virtual real pdf(int state, int action, int next_state) const;
	/// Return the set of next states.
	/// In this case, if a state has not been visited before, then we assume that the next-state set is empty. This means that value iteration will stop upon reaching this state-action pair.
	const DiscreteStateSet& getNextStates(int state, int action) const
	{
		DiscreteStateAction SA(state, action);
		auto got = next_states.find(SA);
		if (got == next_states.end() ){
			return empty_set;
		} else {
			return got->second;
		}
	}

};

#endif
