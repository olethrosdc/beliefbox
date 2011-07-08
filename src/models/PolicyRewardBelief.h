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
#include "Demonstrations.h"
#include "DiscretePolicy.h"
#include "RewardDistribution.h"
#include "RewardBelief.h"
#include "DiscretePolicy.h"
#include "ValueIteration.h"
#include "ExponentialDistribution.h"

#include <vector>


class DiscreteSpaceRewardDistribution;
class ExponentialDistribution;

/** Prior on policy given reward belief */
class PolicyRewardBelief
{
protected:
	DiscreteMDP mdp; ///< the actual MDP (transitions assumed known here)
	int n_states; ///< the number of states
	int n_actions; ///< the number of actions
	ExponentialDistribution softmax_prior; ///< softmax prior
	DirichletRewardBelief reward_prior; ///< reward prior
	real gamma; ///< value of gamma (assumed known here)
    real epsilon; ///< accuracy
	ValueIteration value_iteration;

	// -- the following values are stored during the sampling procedure -- //
	std::vector<Matrix> rewards; ///< set of reward function samples
	std::vector<FixedDiscretePolicy> policies; ///< storage for sampled policies from the belief
	std::vector<real> betas; ///< values of beta sampled
	std::vector<real> sample_counts; ///< amount of times we drew each sample

public:
    PolicyRewardBelief(real lambda,
					   real gamma_,
					   const DiscreteMDP& mdp_);
	virtual ~PolicyRewardBelief();
	
	virtual FixedDiscretePolicy* CalculatePosterior(Demonstrations<int, int>& D);
	virtual real logLikelihood(const Demonstrations<int, int>&D,
							   const  FixedDiscretePolicy& policy) const;

	/// Calculate  $log P(a^T \mid s^T, \pi)$
	virtual real Likelihood(const Demonstrations<int, int>&D,
							const  FixedDiscretePolicy& policy) const
	{
		return exp(logLikelihood(D, policy));
	}

	virtual FixedDiscretePolicy samplePolicy(Matrix& R, real beta);
    virtual FixedDiscretePolicy* getPolicy();
    /// Set accuracy
	 void setAccuracy(real epsilon_)
	{
		epsilon = epsilon_;
		assert(epsilon > 0);
        printf("# setting accuracy to %f\n", 
               epsilon);
	}
	void MHSampler(Demonstrations<int, int>& D, int n_iterations, int n_chains);
	void MonteCarloSampler(Demonstrations<int, int>& D, int n_iterations);
};



#endif
