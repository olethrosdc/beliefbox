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

#include "RewardPolicyBelief.h"
#include <vector>

/// Create from a fixed set of reward distributions
RewardPolicyBelief::RewardPolicyBelief(int n_states, int n_actions,
									   DiscreteMDP& mdp_,
									   const std::vector<DiscreteSpaceRewardDistribution> rewards_)

	: policy_belief(n_states, n_actions),
	  mdp(mdp_),
	  rewards(rewards_)
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
}

/// Enumerate all index reward functions
RewardPolicyBelief::RewardPolicyBelief(int n_states, int n_actions,
									   DiscreteMDP& mdp_)
	: policy_belief(n_states, n_actions),
	  mdp(mdp_)
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			DiscreteSpaceRewardDistribution R_sa(n_states, n_actions);
			R_sa.setFixedReward(s, a, 1.0);
			rewards.push_back(R_sa);
		}
	}
}


/// Calculate a posterior over reward functions
real RewardPolicyBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
	policy_belief.CalculatePosterior(D);
	for (uint i=0; i<rewards.size(); ++i)  {
		
		
	}
	return 1.0;
}


