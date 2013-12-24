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
#include "DiscreteStateSet.h"
#include "StateAction.h"
#include "HashCombine.h"
#include <map>
#include <unordered_map>

template <typename StateType, typename ActionType>
class TransitionDistribution
{
public:
    virtual ~TransitionDistribution() {}
    virtual StateType generate(StateType state, ActionType action) = 0;
    virtual real pdf(StateType state, ActionType action, StateType next_state) = 0;
};

template <typename StateType, typename ActionType>
class Transition
{
public:
	StateType state;
	ActionType action;
	StateType next_state;
	Transition(StateType s, 
			   ActionType a,
			   StateType s2)
		: state(s), action(a), next_state(s2)
	{}
	bool operator==(const Transition& rhs) const
	{
		if ((rhs.state != state)
			|| (rhs.action != action)
			|| (rhs.next_state != next_state)
			) {
			return false;
		}
		return true;
	}
};

typedef Transition<int, int> DiscreteTransition;

namespace std
{
	template<>
 	struct hash<DiscreteTransition>
 	{
		/// The hash is a xor of three hashes
		std::size_t operator() (const DiscreteTransition & s) const
		{
			std::size_t seed = 0;
			hash_combine(seed, hash<int>()(s.state));
			hash_combine(seed, hash<int>()(s.action));
			hash_combine(seed, hash<int>()(s.next_state));
			return seed;
		}
		
 	};
}
 
 
/** Discrete transition distribution.

	In this model, there is no remaining probability mass. So, the probability of going to any particular state is zero.
	
 */
template<>
class TransitionDistribution<int, int>
{
public:
	int n_states; ///< the maximum number of state
	int n_actions; ///< the maximum number of actions
	DiscreteStateSet empty_set; ///< included for convenience
	/// The implementation of the discrete transition distribution
	std::unordered_map<DiscreteTransition, real> P;  ///< gives the actual probabilities
	std::unordered_map<DiscreteStateAction, DiscreteStateSet> next_states; ///< next states for quick access
	TransitionDistribution(int n_states_, int n_actions_)
		: n_states(n_states_),
		  n_actions(n_actions_)
	{
	}
	virtual ~TransitionDistribution() {}

	/// Set a state transition
	virtual void SetTransition(int state, int action, int next_state, real probability);

	/// Get a state transition
	virtual real GetTransition(int state, int action, int next_state) const;

	/// Generate a next state
	virtual int generate(int state, int action) const;
	/// Get the probability of the next state
	virtual real pdf(int state, int action, int next_state) const;
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

typedef TransitionDistribution<int, int> DiscreteTransitionDistribution;
#endif
