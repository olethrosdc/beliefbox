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

// models
#include "RewardPolicyBelief.h"

// algorithms
#include "DiscretePolicy.h"
#include "ValueIteration.h"
#include "PolicyEvaluation.h"

// system
#include <vector>
#include <set>

/// Create from a fixed set of reward distributions
RewardPolicyBelief::RewardPolicyBelief(real lambda_,
                                       real gamma_,
									   DiscreteMDP& mdp_,
									   const std::vector<DiscreteSpaceRewardDistribution> rewards_)

	: n_states(mdp_.getNStates()),
	  n_actions(mdp_.getNActions()),
	lambda(lambda_),
	  policy_belief(n_states, n_actions),
      gamma(gamma_),
	  mdp(mdp_),
	  rewards(rewards_),
	  P_rewards(rewards.size())
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
    assert(gamma >= 0 && gamma <= 1);
    setAccuracy(1e-3);
}

/// Enumerate all index reward functions
RewardPolicyBelief::RewardPolicyBelief(real lambda_,
                                       real gamma_,
									   DiscreteMDP& mdp_)
	: n_states(mdp_.getNStates()),
	  n_actions(mdp_.getNActions()),
	  lambda(lambda_),
      policy_belief(n_states, n_actions),
      gamma(gamma_),
	  mdp(mdp_),
	  P_rewards(rewards.size())
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
    assert(gamma >= 0 && gamma <= 1);
    setAccuracy(1e-3);
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
	real epsilon = 1e-3; ///< minimu precision
	int max_iter = 10e4; ///< maximum number of iterations

	policy_belief.CalculatePosterior(D);
	//--  resample from the belief -- //
	
	// reset policies vector
	for (uint i=0; i<policies.size(); ++i) {
		delete policies[i];
	}
	policies.resize(n_samples);

	// add new samples
	for (int i=0; i<n_samples; ++i) {
		policies[i] = policy_belief.Sample();
	}


	//-- calculate probability of each policy --//
	int n_rewards = P_rewards.Size();
	assert(n_rewards == rewards.size());

	//-- Make a matrix that contains the optimality of each policy --//
	P_rewards.Clear();
	Matrix L(n_samples, n_rewards);
	std::set<real> loss_vector;
	for (int i=0; i<P_rewards.Size(); ++i) {
		// Change MDP reward
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                mdp.setFixedReward(s, a, rewards[i].expected(s, a));
            }
        }		
		// Calculate value of optimal policy
		ValueIteration VI(&mdp, gamma);
		VI.ComputeStateValues(epsilon, max_iter);
		for (int j=0; j<n_samples; ++j) {
			// Calculate value of actual policy;
			PolicyEvaluation PE(policies[j], &mdp, gamma);
			PE.ComputeStateValues(epsilon);
			L(i, j) = VI.getValue(0) - PE.getValue(0);
			for (int s=1; s<n_states; ++s) {
				real DV_s = VI.getValue(s) - PE.getValue(s);
				if (DV_s > L(i, j)) {
					L(i, j) = DV_s;
				}
			}
		}
	}
	
	return 1.0;
}


/// Virtual destructor
RewardPolicyBelief::~RewardPolicyBelief()
{
	for (uint i=0; i<policies.size(); ++i) {
		delete policies[i];
	}	
}

