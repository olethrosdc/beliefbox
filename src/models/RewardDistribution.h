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

#include "real.h"

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
    real default_reward;
    int goal_state;
    real goal_reward;
public:
    DiscreteSpaceRewardDistribution();
    virtual ~DiscreteSpaceRewardDistribution();
    virtual real generate(int state, int action) const;
    virtual real expected(int state, int action) const;
    virtual real pdf(int state, int action, real reward) const;
};



#endif
