// -*- Mode: c++ -*-
// copyright (c) 2005-2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_MDP_H
#define DISCRETE_MDP_H

#ifndef NDEBUG
#include "SmartAssert.h"
#endif
#include "MDP.h"
#include "Distribution.h"
#include "RewardDistribution.h"
#include "TransitionDistribution.h"
#include "Matrix.h"
#include "DiscreteStateSet.h"
#include <vector>
#include <set>




template<>
class MDP<int, int> {
protected:
	real reward; ///< current reward
	int state; ///< current state;
	int n_states; ///< number of states (or dimensionality of state space)
	int n_actions; ///< number of actions (or dimensionality of action space)
	int N;
	inline int getID (int s, int a) const
	{
#ifndef NDEBUG
		//DISABLED_ASSERT(s>=0 && s<n_states)(s)(n_states);
		//DISABLED_ASSERT(a>=0 && a<n_actions)(a)(n_actions);
#endif
		return s*n_actions + a;
	}

public:
	/// Reward distribution
	DiscreteSpaceRewardDistribution reward_distribution;
	/// Transition distribution
	DiscreteTransitionDistribution transition_distribution; 

	/// Default constructor
	MDP<int, int>(int n_states_, int n_actions_,
				  real** initial_transitions = NULL);
	/// Copy constructor
	MDP<int,int> (const MDP<int,int>& mdp);
	/// Mean MDP constructor
	MDP<int,int> (const std::vector<const MDP<int,int>*> &mdp_list,
				  const Vector& w);

	inline int StateLowerBound() const
	{
		return 0;
	}
	inline int StateUpperBound() const
	{
		return n_states;
	}
	inline int getNStates() const
	{
		return n_states;
	}
	inline int getNActions() const
	{
		return n_actions;
	}
	virtual ~MDP<int, int>();
	int getState() const
	{
		return state;
	}
	real getReward() const
	{
		return reward;
	}
	int Reset(int new_state)
	{
		reward = 0.0;
		return setState(new_state);
	}
	int setState(int new_state)
	{
		assert(new_state >=0 && new_state < n_states);
		state = new_state;
		return state;
	}
	// generate a new state given the current state and action, then set the current state to be the new state.
	real Act (int a)
	{
		reward = generateReward(state, a);
		state = generateState(state, a);
		return reward;
	}
	/// Simply show the model
	virtual void ShowModel() const;

	/// Output the model in a dot file for visualizaton
	virtual void dotModel(FILE* fout) const;

	/// Generate a reward from the model
	real generateReward (const int& s, const int& a) const;

	/// Generate a next state from the model
	int generateState (const int& s, const int& a) const;

	/// Get the reward probability 
	virtual real getRewardProbability (const int& s, const int& a, real r) const
	{
		return reward_distribution.pdf(s, a, r);
	}

	virtual real getTransitionProbability (const int& s, const int& a, const int& s2) const
	{
		//int ID = getID (s, a);                
		//assert (s2>=0 && s2<n_states);
		//return P(ID, s2);
		return transition_distribution.pdf(s, a, s2);
	}
	virtual real getExpectedReward (const int& s, const int& a) const
	{
		return reward_distribution.expected(s,a);
	}
	virtual void setTransitionProbability(int s, int a, int s2, real p)
	{
		transition_distribution.SetTransition(s, a, s2, p);
	}
	virtual void setTransitionProbabilities(int s, int a, const Vector& p, real threshold = 0)
	{
		assert(s>=0 && s<n_states);
		assert(p.Size() == n_states);
		for (int s2=0; s2<n_states; ++s2) {
			transition_distribution.SetTransition(s, a, s2, p(s2));
		}
	}
	virtual const DiscreteStateSet& getNextStates(int s, int a) const
	{
		return transition_distribution.getNextStates(s, a);
	}

	void AperiodicityTransform(real tau);
	bool Check() const;
	real CalculateDiameter() const;
	virtual void setRewardDistribution(int s, int a, Distribution* reward)
	{
		assert(s>=0 && s<n_states);
		reward_distribution.setRewardDistribution(s, a, reward);
	}
	virtual void addRewardDistribution(int s, int a, Distribution* reward)
	{
		assert(s>=0 && s<n_states);
		reward_distribution.addRewardDistribution(s, a, reward);
	}
	virtual void addFixedReward(int s, int a, real reward)
	{
		assert(s>=0 && s<n_states);
		reward_distribution.addFixedReward(s, a, reward);
	}
	virtual void setFixedReward(int s, int a, real reward)
	{
		assert(s>=0 && s<n_states);
		reward_distribution.setFixedReward(s, a, reward);
	}
	virtual void setFixedRewards(const Matrix& R)
	{
		assert(R.Rows() == n_states);
		assert(R.Columns() == n_actions);
		for (int s=0; s<n_states; ++s) {
			for (int a=0; a<n_actions; ++a) {
				setFixedReward(s, a, R(s,a));
			}
		}
	}
};

typedef MDP<int, int> DiscreteMDP;

#endif
