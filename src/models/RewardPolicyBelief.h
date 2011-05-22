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
    real gamma;
    DirichletProductPolicyBelief policy_belief;
	DiscreteMDP mdp;
	std::vector<DiscreteSpaceRewardDistribution> rewards;
    int n_samples;
public:
    RewardPolicyBelief(int n_states, int n_actions,
                       real gamma_,
					   DiscreteMDP& mdp_,
					   const std::vector<DiscreteSpaceRewardDistribution> rewards_);	

    RewardPolicyBelief(int n_states, int n_actions,
                       real gamma_,
					   DiscreteMDP& mdp_);

	/// Virtual destructor
	virtual ~RewardPolicyBelief()
	{ }
	
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
