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
									   const DiscreteMDP& mdp_,
									   const std::vector<DiscreteSpaceRewardDistribution*> rewards_)

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
    setAccuracy(1e-1);
}

/// Enumerate all index reward functions
RewardPolicyBelief::RewardPolicyBelief(real lambda_,
                                       real gamma_,
									   const DiscreteMDP& mdp_)
	: n_states(mdp_.getNStates()),
	  n_actions(mdp_.getNActions()),
	  lambda(lambda_),
      policy_belief(n_states, n_actions),
      gamma(gamma_),
	  mdp(mdp_)
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
    assert(gamma >= 0 && gamma <= 1);
    setAccuracy(1e-1);

	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			DiscreteSpaceRewardDistribution* R_sa = new DiscreteSpaceRewardDistribution(n_states, n_actions);
            R_sa->setFixedReward(s, a, 1.0); 
			rewards.push_back(R_sa);
		}
	}
    P_rewards.Resize(rewards.size());
}

/// Sample a set of reward functions
RewardPolicyBelief::RewardPolicyBelief(real lambda_,
                                       real gamma_,
									   const DiscreteMDP& mdp_,
                                       const DirichletDistribution& reward_prior,
                                       int n_reward_samples)
	: n_states(mdp_.getNStates()),
	  n_actions(mdp_.getNActions()),
	  lambda(lambda_),
      policy_belief(n_states, n_actions),
      gamma(gamma_),
	  mdp(mdp_)
{
	assert(n_states == mdp.getNStates());
	assert(n_actions == mdp.getNActions());
    assert(gamma >= 0 && gamma <= 1);
    setAccuracy(1e-1);

    for (int i=0; i<n_reward_samples; ++i) {
        DiscreteSpaceRewardDistribution* R_sa = new DiscreteSpaceRewardDistribution(n_states, n_actions);
        Vector R = reward_prior.generate();
        int i = 0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a, ++i) {
                R_sa->setFixedReward(s, a, R(i)); 
            }
        }
        rewards.push_back(R_sa);
    }
    P_rewards.Resize(rewards.size());
}



/// Calculate a posterior over reward functions
DiscretePolicy* RewardPolicyBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
	real epsilon = 1e-3; ///< minimum precision
	int max_iter = 1e3; ///< maximum number of iterations

	// --------  resample from the belief -------- //
	// reset policies vector
	for (uint i=0; i<policies.size(); ++i) {
		delete policies[i];
	}
	policies.resize(n_policies);

	// add new samples
	policy_belief.CalculatePosterior(D);
	for (int i=0; i<n_policies; ++i) {
		policies[i] = policy_belief.Sample();
	}


	// -------- calculate probability of each policy -------- //
	int n_rewards = P_rewards.Size();
	assert(n_rewards == (int) rewards.size());

	// -------- create loss matrix -------- //
	P_rewards.Clear();
	Matrix L(n_rewards, n_policies);
	std::set<real> loss_vector;

	printf ("# calculating %d x %d loss matrix\n", n_rewards, n_policies);
	for (int i=0; i<P_rewards.Size(); ++i) {
		//mdp.Check();
		//mdp.ShowModel();
		// Change MDP reward to the i-th reward
        #if 1
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                mdp.setFixedReward(s, a, rewards[i]->expected(s, a));
				/// BUG - for some reason the MDP value iteration does not work
                //mdp.setFixedReward(s, a, urandom());
				//printf ("R(%d, %d) = %f (%f)\n", s, a, rewards[i]->expected(s, a),
				//mdp.getExpectedReward(s, a));
            }
        }	
        #else
        mdp.reward_distribution = *rewards[i];
        #endif
		//mdp.ShowModel();
		// Calculate value of optimal policy for the i-th reward function
		ValueIteration VI(&mdp, gamma);
		VI.ComputeStateActionValues(epsilon, max_iter);

		// Calculate the loss for each policy sample
		for (int j=0; j<n_policies; ++j) {
			// Calculate value of actual policy;
			PolicyEvaluation PE(policies[j], &mdp, gamma);
			PE.ComputeStateValues(epsilon);
			L(i, j) = VI.getValue(0) - PE.getValue(0);
			for (int s=0; s<n_states; ++s) {
				real DV_s = VI.getValue(s) - PE.getValue(s);
				//printf ("# s: %d, V(s)=%f, Vk(s)=%f\n", s, VI.getValue(s), PE.getValue(s));
				if (DV_s > L(i, j)) {
					L(i, j) = DV_s;
				}
			}
			//printf ("Inserting L(%d, %d) = %f\n", i, j, L(i, j));
			loss_vector.insert(L(i,j));
		}
	}

	
	// -------- calculate the final posterior -------- //
	int n_losses = loss_vector.size();
	Vector loss_prior (n_losses);

	// first calculate the loss measures
	{
        printf("# loss measure\n");
		std::set<real>::iterator it;
		int k = 0;
		real prev_loss = 0.0;
		for (it=loss_vector.begin(); it!=loss_vector.end(); ++it, ++k) {
			real loss = *it;
			loss_prior(k) = exp(-lambda * prev_loss) - exp(-lambda * loss);
			prev_loss = loss;
		}
	}

	// now calculate the posterior measure over policies
    printf("# reward posterior (%d x %d x %d = %d)\n",
           n_rewards,
           n_policies,
           loss_vector.size(),
           n_rewards * n_policies * loss_vector.size());
	Vector reward_posterior(n_rewards);
	reward_posterior.Clear();
	for (int i=0; i<n_rewards; ++i) {
		for (int j=0; j<n_policies; ++j) {
            #if 0
			std::set<real>::iterator it;
			int k = 0;
			for (it=loss_vector.begin(); it!=loss_vector.end(); ++it, ++k) {
				real loss = *it;
				if (L(i, j) <= loss) {
					reward_posterior(i) += loss_prior(k);
				}
			}
            #else
            reward_posterior(i) += exp(-lambda * L(i,j));
			#endif
		}
	}

    reward_posterior /= reward_posterior.Sum();
    //reward_posterior.print(stdout); printf ("# Reward posterior!\n");
    int arg_max = ArgMax(reward_posterior);
    rewards[arg_max]->getExpectedRewardVector().print(stdout); printf(" # MAP\n");
    Vector R(n_states * n_actions);
    R.Clear();
    for (int i=0; i<n_rewards; ++i) {
        R += rewards[i]->getExpectedRewardVector() * reward_posterior(i);
    }
    R.print(stdout); printf("# Expected\n");
    {
        int k=0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a, ++k) {
                mdp.setFixedReward(s, a, R(k));
            }
        }	
        ValueIteration VI(&mdp, gamma);
		VI.ComputeStateActionValues(epsilon, max_iter);
        return VI.getPolicy();
    }
}


/// Virtual destructor
RewardPolicyBelief::~RewardPolicyBelief()
{
	for (uint i=0; i<policies.size(); ++i) {
		delete policies[i];
	}	
}

