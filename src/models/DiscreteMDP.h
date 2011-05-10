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

#include "SmartAssert.h"
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
    int state;
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

    inline int GetNStates() const
    {
        return n_states;
    }
    inline int GetNActions() const
    {
        return n_actions;
    }
    virtual ~MDP<int, int>();
	int getState() const
	{
		return state;
	}
    // generate a new state given the current state and action, then set the current state to be the new state.
    real Act (int a)
    {
        real r = generateReward(state, a);
        state = generateState(state, a);
        return r;
    }
    virtual void ShowModel() const;
    virtual void dotModel(FILE* fout) const;
    real generateReward (int s, int a) const;
    int generateState (int s, int a) const;
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
		reward_distribution.setRewardDistribution(s, a, reward);
	}
	inline void addRewardDistribution(int s, int a, Distribution* reward)
	{
		reward_distribution.addRewardDistribution(s, a, reward);
	}
	inline void addFixedReward(int s, int a, real reward)
	{
		reward_distribution.addFixedReward(s, a, reward);
	}
	inline void setFixedReward(int s, int a, real reward)
	{
		reward_distribution.setFixedReward(s, a, reward);
	}
};

typedef MDP<int, int> DiscreteMDP;

#endif
