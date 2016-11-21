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

class DiscreteJointAction
{
public:
	int n_players;
	int n_actions;
	std::vector<int> action;
	DiscreteJointAction(int n_players_,
						int n_actions_)
		: n_players(n_players_),
		  n_actions(n_actions_),
		  action(n_players_)
	{
		for (int i=0; i<n_players; ++i) {
			action[i] = 0;
		}
	}
	
	/// Convert the action vector to an integer
	int to_int()
	{
		int F = 1;
		int a = 0;
		for (int k=0; k<n_players; ++k) {
			a += action[k];
			F *= n_actions;
		}
		return a;
	}
};

// maybe use a vector of MoveType instead of that

class DiscreteMultiTransitionDistribution
{
public:
	int n_players;
	int n_actions;
	int n_joint_actions;
	DiscreteTransitionDistribution dist;
	DiscreteTransitionDistribution(n_players_, n_actions_)
		: n_players(n_players_),
		  n_actions(n_actions_)
	{
		n_joint_actions = n_actions;
		for (int k=1; k<n_players; ++k) {
			n_joint_actions *= n_actions;
		}
		dist = DiscreteTransitionDistribution(n_states, n_joint_actions);
	}
	real pdf(int s, DiscreteJointAction& a, int s2)
	{
		return dist.pdf(s, a.to_int(), s2);
	}
};

// This is just a simple template for Multi-Agent MDPs
template<>
class MMDP<int, int> {
{
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
public:
    DiscreteMMDP(int n_players_,
				 int n_states_,
				 int n_actions_)
		: n_players(n_players_),
		  n_actions(n_actions_),
		n_states(n_states_),
		state(0),
		action(n_players),
		transition_distribution(n_states, n_players, n_actions),
		simultaneous_moves(false),
		reward(n_states)
		{
			// nothing to do in the constructor itself
			
		}
	virtual ~MMDP() {}

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
    virtual real getReward (const int& s,
							const DiscreteJointAction& a) const
    {
        return reward[s];
    }

	/// Generate a state from the transition distribution but don't change the MMDP state.x
    virtual StateType generateState(const int& s,
									const DiscreteJointAction& a) const
    {
        return transition_distribution.generate(s, a);
    }
    /// generate a new state given the current state and action, then set the current state to be the new state.
    virtual real Act (DiscreteJointAction& a)
    {
        real r = generateReward(state, a);
        state = generateState(state, a);
        return r;
    }
	/// Generate a state from the transition distribution for the current state but don't change the MMDP state.
    virtual StateType generateState(const ActionType& a)
    {
        return generateState(state, a);
    }
};

/// Short-hand for this class
typedef MMDP<int, int> DiscreteMMDP;
	
#endif
	
