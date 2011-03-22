// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "RewardDistribution.h"

// -- Goal state reward distribution -- //

/// Constructor. Set up the default reward to be received
real GoalStateRewardDistribution(real default_reward_) :
    default_reward(default_reward_),
    goal_state(goal_state_),
    goal_reward(goal_reward_)
{
    // empty
}

/// Destructor - nothing to do
GoalStateRewardDistribution::~GoalStateRewardDistribution() 
{
}

/// Generate a value
real GoalStateRewardDistribution::generate(int state, int action)
{
    if (state == goal_state) {
        return goal_reward;
    }
    return default_reward;
}

real GoalStateRewardDistribution::expected(int state, int action)
{
    if (state == goal_state) {
        return goal_reward;
    }
    return default_reward;
}

real GoalStateRewardDistribution::pdf(int state, int action, real reward)
{
    real r = generate(state, action);
    if (fabs(r - reward) < 1e-6) {
        return 1.0;
    } else {
        return 0.0;
    }
    
}
