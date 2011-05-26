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
#include "Distribution.h"

// -- Goal state reward distribution -- //

/// Constructor. Set up the default reward to be received
GoalStateRewardDistribution::GoalStateRewardDistribution(real default_reward_, 
														 int goal_state_,
														 real goal_reward_)
	: default_reward(default_reward_),
	  goal_state(goal_state_),
	  goal_reward(goal_reward_)
{
    // empty
}

/// Destructor - nothing to do
GoalStateRewardDistribution::~GoalStateRewardDistribution() 
{
}

/// Generate a value - always the expected value
real GoalStateRewardDistribution::generate(int state, int action) const
{
	return expected(state, action);
}

/// Get expected value
real GoalStateRewardDistribution::expected(int state, int action) const
{
    if (state == goal_state) {
        return goal_reward;
    }
    return default_reward;
}

/// Generate pdf -- in fact, a probability since there are just two outcomes.
real GoalStateRewardDistribution::pdf(int state, int action, real reward) const
{
    real r = expected(state, action);
    if (fabs(r - reward) < 1e-6) {
        return 1.0;
    } else {
        return 0.0;
    }
    
}

// -- General reward distribution for discrete environments -- //

/// Constructor. Set up the default reward to be received
DiscreteSpaceRewardDistribution::DiscreteSpaceRewardDistribution(int n_states_, int n_actions_)
    : n_states(n_states_),
      n_actions(n_actions_),
	  R(n_states * n_actions),
	  ER(n_states * n_actions) 
{
    // empty
	for (uint i=0; i<R.size(); ++i) {
		R[i] = NULL;
	}
	ER.Clear();
    //printf ("%d %d %d %d\n", n_states, n_actions, R.size(), ER.Size());
}

/// Copy constructor. Do not actually copy anything!
DiscreteSpaceRewardDistribution::DiscreteSpaceRewardDistribution(const DiscreteSpaceRewardDistribution& rhs)
{
	n_states = rhs.n_states;
	n_actions = rhs.n_actions;
	ER = rhs.ER;
}

/// Assignment operator. Do not copy anything!
DiscreteSpaceRewardDistribution& DiscreteSpaceRewardDistribution::operator= (const DiscreteSpaceRewardDistribution& rhs)
{
	if (this == &rhs) return *this;
	n_states = rhs.n_states;
	n_actions = rhs.n_actions;
	ER = rhs.ER;
	return *this;
}



/// Destructor - remove all distribution vectors
DiscreteSpaceRewardDistribution::~DiscreteSpaceRewardDistribution() 
{
    // clear all distributions that have been added
    for (uint i = 0; i<distribution_vector.size(); ++i) {
        if (distribution_vector[i]) {
            delete distribution_vector[i];
			distribution_vector[i] = NULL;
        }
    }

}

/// Generate a value - always the expected value
real DiscreteSpaceRewardDistribution::generate(int state, int action) const
{
	int ID = getID(state, action);
	if (R[ID]) {
		return R[ID]->generate();
	} else {
		return ER[ID];
	}
}

/// Get expected value
real DiscreteSpaceRewardDistribution::expected(int state, int action) const
{
	int ID = getID(state, action);
	if (R[ID]) {
		return R[ID]->getMean();
	} else {
		return ER[ID];
	}
}

/// Generate pdf.
real DiscreteSpaceRewardDistribution::pdf(int state, int action, real reward) const
{
	int ID = getID (state, action);
	if (R[ID]) {
		return R[ID]->pdf(reward);
	}
    real r = expected(state, action);
    if (fabs(r - reward) < 1e-6) {
        return 1.0;
    } else {
        return 0.0;
    }
    
}


void DiscreteSpaceRewardDistribution::setRewardDistribution(int s, int a, Distribution* reward)
{   
	int ID = getID (s, a);
	R[ID] = reward;
	ER(ID) = reward->getMean();
}

// only use this function once per state-action pair
void DiscreteSpaceRewardDistribution::addRewardDistribution(int s, int a, Distribution* reward)
{   
	int ID = getID (s, a);
    assert(reward);
	distribution_vector.push_back(reward);
	R[ID] = reward;
	ER(ID) = reward->getMean();
}
// only use this function once per state-action pair
void DiscreteSpaceRewardDistribution::addFixedReward(int s, int a, real reward)
{   
	SingularDistribution* distribution = new SingularDistribution(reward);
	addRewardDistribution(s, a, distribution);
}

// you can use this function more than once per state-action pair
void DiscreteSpaceRewardDistribution::setFixedReward(int s, int a, real reward)
{   
	int ID = getID (s, a);

	if (R[ID]) {
		R[ID]->setMean(reward);
		ER[ID] = reward;
	} else {
		SingularDistribution* distribution = new SingularDistribution(reward);
		addRewardDistribution(s, a, distribution);
		ER[ID] = reward;
	}
	
}
