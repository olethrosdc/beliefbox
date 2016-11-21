// -*- Mode: c++ -*-
// copyright (c) 2016 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_MMDP_H
#define DISCRETE_MMDP_H


#include <vector>
#include <memory>
#include <cassert>

#include "MMDP.h"
#include "TransitionDistribution.h"
#include "JointAction.h"
#include "Random.h"

// maybe use a vector of MoveType instead of that

class DiscreteMultiTransitionDistribution
{
public:
	int n_players;
	int n_states;
	int n_actions;
	int n_joint_actions;
	std::shared_ptr<DiscreteTransitionDistribution> dist;
	DiscreteMultiTransitionDistribution(const int& n_players_,
										const int& n_states_,
										const int& n_actions_)
		: n_players(n_players_),
		  n_states(n_states_),
		  n_actions(n_actions_)
	{
		n_joint_actions = n_actions;
		for (int k=1; k<n_players; ++k) {
			n_joint_actions *= n_actions;
		}
		dist = std::make_shared<DiscreteTransitionDistribution>(n_states, n_joint_actions);
	}
	real pdf(const int& s, const DiscreteJointAction& a, const int& s2) const
	{
		return dist->pdf(s, a.to_int(), s2);
	}
	int generate(const int& s, const DiscreteJointAction& a) const
	{
		return dist->generate(s, a.to_int());
	}
	void setTransitionProbability(const int& s, const DiscreteJointAction& a, const int& s2, real p)
	{
		dist->SetTransition(s, a.to_int(), s2, p);
	}
};

// This is just a simple template for Multi-Agent MDPs
template<>
class MMDP<int, int> {
protected:
	int n_players;
	int n_actions;
	int n_states;
    int state;
	DiscreteJointAction action;
	DiscreteMultiTransitionDistribution transition_distribution;
	int current_player; ///< We need to know who is playing
	bool simultaneous_moves; ///< should the players move simultaneously
	std::vector<real> reward; ///< a simple reward vector for now
	int n_joint_actions;
public:
    MMDP<int,int>(int n_players_,
				 int n_states_,
				 int n_actions_)
		: n_players(n_players_),
		  n_actions(n_actions_),
		n_states(n_states_),
		state(0),
		action(n_players, n_actions),
		transition_distribution(n_states, n_players, n_actions),
		simultaneous_moves(false),
		reward(n_states)
		{
			// Calculate the number of joint actions
			n_joint_actions = n_actions;
			for (int k=1; k<n_players; ++k) {
				n_joint_actions *= n_actions;
			}
		}
	virtual ~MMDP<int,int>() {}

	int getNJointActions()
	{
		return n_joint_actions;
	}
	/// Get the current stat
	int getState() const
	{
		return state;
	}
	/// Get the transition probability
    virtual real getTransitionProbability (const int& s,
                                           const DiscreteJointAction& a,
                                           const int& s2) const
    {
        return transition_distribution.pdf(s, a, s2);
    }
	virtual void setTransitionProbability (const int& s,
                                           const DiscreteJointAction& a,
                                           const int& s2,
										   real p) 
    {
        transition_distribution.setTransitionProbability(s, a, s2, p);
    }

    virtual real getReward (const int& s) const
    {
        return reward[s];
    }
	virtual void setReward (const int& s, const real& r) 
    {
        reward[s] = r;
    }

	/// Generate a state from the transition distribution but don't change the MMDP state.x
    virtual int generateState(const int& s,
									const DiscreteJointAction& a) const
    {
        return transition_distribution.generate(s, a);
    }
    /// generate a new state given the current state and action, then set the current state to be the new state.
    virtual real Act (DiscreteJointAction& a)
    {
        real r = getReward(state);
        state = generateState(state, a);
        return r;
    }
	/// Generate a state from the transition distribution for the current state but don't change the MMDP state.
    virtual int generateState(const DiscreteJointAction& a)
    {
        return generateState(state, a);
    }
};

/// Short-hand for this class
typedef MMDP<int, int> DiscreteMMDP;
	
#endif
	
