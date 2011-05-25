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

#ifndef REWARD_POLICY_BELIEF_H
#define REWARD_POLICY_BELIEF_H

#include "PolicyBelief.h"
#include "DiscreteMDP.h"

/** Prior on epsilon optimality of policies */
class RewardPolicyBelief
{
protected:
	real lambda; ///< Exponential distribution parameter for the sub-optimality of policies
    DirichletProductPolicyBelief policy_belief;   ///< prior about policies
	real gamma; ///< value of gamma (assumed known here)
	DiscreteMDP mdp; ///< the actual MDP (transitions assumed known here)
	std::vector<DiscreteSpaceRewardDistribution> rewards; ///< set of reward functions
    int n_samples; ///< number of samples required
	std::vector<DiscretePolicy*> policies; ///< storage for sampled policies from the belief
	Vector P_rewards; ///< posterior probability of each reward function
public:
    RewardPolicyBelief(int n_states, int n_actions,
					   real lambda_,
                       real gamma_,
					   DiscreteMDP& mdp_,
					   const std::vector<DiscreteSpaceRewardDistribution> rewards_);	

    RewardPolicyBelief(int n_states, int n_actions,
					   real lambda_,
                       real gamma_,
					   DiscreteMDP& mdp_);

	virtual ~RewardPolicyBelief();
	
	virtual real CalculatePosterior(Demonstrations<int, int>& D);
	
	/// Set number of samples
	void setNSamples(int n_samples_)
	{
		n_samples = n_samples_;
	}

	/// Set accuracy
	void setAccuracy(real epsilon)
	{
        n_samples = (int) ceil(pow((1 - gamma) * epsilon, -2.0));
        printf("# setting accuracy to %f -> n_samples = %d\n", 
               epsilon,
               n_samples);
	}
};

#endif
