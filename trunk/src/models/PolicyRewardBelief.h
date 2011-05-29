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

#ifndef POLICY_REWARD_BELIEF_H
#define POLICY_REWARD_BELIEF_H

#include "DiscreteMDP.h"

/** Prior on policy given belief */
class RewardPolicyBelief
{
protected:
	int n_states; ///< the number of states
	int n_actions; ///< the number of actions
	real lambda; ///< Exponential distribution parameter for the temperature
	real gamma; ///< value of gamma (assumed known here)
    real epsilon; ///< accuracy
	DiscreteMDP mdp; ///< the actual MDP (transitions assumed known here)
	std::vector<DiscreteSpaceRewardDistribution*> rewards; ///< set of reward function samples
	std::vector<DiscretePolicy*> policies; ///< storage for sampled policies from the belief
	Vector P_rewards; ///< posterior probability of each reward function
    Vector estimated_reward;
public:
    PolicyReward Belief(real lambda_,
                        real gamma_,
                        const DiscreteMDP& mdp_,
                        const std::vector<DiscreteSpaceRewardDistribution*> rewards_);	
	virtual ~RewardPolicyBelief();
	
	virtual DiscretePolicy* CalculatePosterior(Demonstrations<int, int>& D);
	
	/// Set accuracy
	void setAccuracy(real epsilon_)
	{
		epsilon = epsilon_;
		assert(epsilon > 0);
        printf("# setting accuracy to %f\n", 
               epsilon);
	}
};



#endif
