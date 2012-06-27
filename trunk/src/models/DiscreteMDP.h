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
#include "Matrix.h"
#include <vector>
#include <set>

typedef std::set<int> DiscreteStateSet;
typedef std::set<int>& DiscreteStateSetRef;




template<>
class MDP<int, int> {
protected:
  real reward; ///< current reward
    int state; ///< current state;
    int n_states; ///< number of states (or dimensionality of state space)
    int n_actions; ///< number of actions (or dimensionality of action space)
    Matrix P; ///< transition distribution
    std::vector<DiscreteStateSet> next_states;
    int N;
    inline int getID (int s, int a) const
    {
#ifndef NDEBUG
        SMART_ASSERT(s>=0 && s<n_states)(s)(n_states);
        SMART_ASSERT(a>=0 && a<n_actions)(a)(n_actions);
#endif
        return s*n_actions + a;
    }

public:
    DiscreteSpaceRewardDistribution reward_distribution;
    
    MDP<int, int>(int n_states_, int n_actions_,
                  real** initial_transitions = NULL);
	//Distribution** initial_rewards = NULL);
    MDP<int,int> (const MDP<int,int>& mdp);
    MDP<int,int> (const std::vector<const MDP<int,int>*> &mdp_list,
                  const Vector& w);

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
    virtual void ShowModel() const;
    virtual void dotModel(FILE* fout) const;
    real generateReward (int s, int a) const;
    int generateState (int s, int a) const;
    inline real getRewardProbability (int s, int a, real r) const
    {
        return reward_distribution.pdf(s, a, r);
    }

    inline real getTransitionProbability (int s, int a, int s2) const
    {
        int ID = getID (s, a);                
        assert (s2>=0 && s2<n_states);
        return P(ID, s2);
    }
    inline real getExpectedReward (int s, int a) const
    {
        return reward_distribution.expected(s,a);
    }
    inline void setTransitionProbability(int s, int a, int s2, real p)
    {
        assert(s>=0 && s<n_states);
        int ID = getID (s, a);
        assert(s2>=0 && s2<n_states);
		P(ID, s2) = p;
        DiscreteStateSet& next = next_states[ID];
        if (p==0) {
            next.erase(s2);
        } else {
            next.insert(s2);
        }
    }
    inline void setTransitionProbabilities(int s, int a, const Vector& p)
    {
        assert(s>=0 && s<n_states);
        int ID = getID (s, a);
        assert(p.Size() == n_states);
        for (int s2=0; s2<n_states; ++s2) {
            P(ID, s2) = p(s2);
            DiscreteStateSet& next = next_states[ID];
            if (p(s2)==0) {
                next.erase(s2);
            } else {
                next.insert(s2);
            }
        }
    }
    inline const DiscreteStateSet& getNextStates(int s, int a) const
    {
        int ID = getID (s,a);
        return next_states[ID];
    }

    void AperiodicityTransform(real tau);
    bool Check() const;
    real CalculateDiameter() const;
	inline void setRewardDistribution(int s, int a, Distribution* reward)
	{
        assert(s>=0 && s<n_states);
		reward_distribution.setRewardDistribution(s, a, reward);
	}
	inline void addRewardDistribution(int s, int a, Distribution* reward)
	{
        assert(s>=0 && s<n_states);
		reward_distribution.addRewardDistribution(s, a, reward);
	}
	inline void addFixedReward(int s, int a, real reward)
	{
        assert(s>=0 && s<n_states);
		reward_distribution.addFixedReward(s, a, reward);
	}
	inline void setFixedReward(int s, int a, real reward)
	{
        assert(s>=0 && s<n_states);
		reward_distribution.setFixedReward(s, a, reward);
	}
	inline void setFixedRewards(const Matrix& R)
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
