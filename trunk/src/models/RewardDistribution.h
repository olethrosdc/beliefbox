// -*- Mode: c++ -*-
// copyright (c) 2007-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REWARD_DISTRIBUTION_H
#define REWARD_DISTRIBUTION_H

#include "Vector.h"
#include "real.h"
#include <vector>

class Distribution;


/** Generic reward distribution */
template <typename StateType, typename ActionType>
class RewardDistribution
{
public:
    /// Empty destructor
    virtual ~RewardDistribution() {}
    /// generate reward
    virtual real generate(StateType state, ActionType action) const = 0;
    /// get expected value
    virtual real expected(StateType state, ActionType action) const = 0;
    /// get pdf
    virtual real pdf(StateType state, ActionType action, real reward) const = 0;
};

/** Goal state reward distribution */
class GoalStateRewardDistribution : public RewardDistribution<int, int>
{
protected:
    real default_reward;
    int goal_state;
    real goal_reward;
public:
    GoalStateRewardDistribution(real default_reward_, int goal_state_, real goal_reward_);
    
    virtual ~GoalStateRewardDistribution();
    virtual real generate(int state, int action) const;
    virtual real expected(int state, int action) const;
    virtual real pdf(int state, int action, real reward) const;
};

/** Generate reward distribution for a discrete space */
class DiscreteSpaceRewardDistribution : public RewardDistribution<int, int>
{
protected:
	int n_states; ///< number of states
	int n_actions; ///< number of actions
    std::vector<Distribution*> R; ///< reward distribution
	std::vector<Distribution*> distribution_vector; ///< for malloc
	Vector ER; ///< expected reward
    inline int getID (int s, int a) const
    {
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }

public:
    DiscreteSpaceRewardDistribution(int n_states_, int n_actions_);
    DiscreteSpaceRewardDistribution(const DiscreteSpaceRewardDistribution& rhs);
    DiscreteSpaceRewardDistribution& operator= (const DiscreteSpaceRewardDistribution& rhs);
    virtual ~DiscreteSpaceRewardDistribution();
    virtual real generate(int state, int action) const;
    virtual real expected(int state, int action) const;
    virtual real pdf(int state, int action, real reward) const;
	void setRewardDistribution(int s, int a, Distribution* reward);
	void addRewardDistribution(int s, int a, Distribution* reward);
	void addFixedReward(int s, int a, real reward);
	void setFixedReward(int s, int a, real reward);
    void Show();
    Vector getExpectedRewardVector() const
    {
        return ER;
    }
};





#endif
