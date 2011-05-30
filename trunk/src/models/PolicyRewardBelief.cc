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

#include "PolicyRewardBelief.h"



/// Create teh blief
PolicyRewardBelief::PolicyRewardBelief(real lambda,
									   real gamma_,
									   const DiscreteMDP& mdp_)
	: mdp(mdp_),
	  n_states(mdp.getNStates()),
	  n_actions(mdp.getNActions()),
	  softmax_prior(lambda),
	  reward_prior(n_states, n_actions),
	  gamma(gamma_),
	  value_iteration(&mdp, gamma)
{
	setAccuracy(1e-3);
}

/// Destroy the belief
PolicyRewardBelief::~PolicyRewardBelief()
{
}

/// This samples the posterior belief given the original sample....
DiscretePolicy* PolicyRewardBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
	
	return NULL;
}


/// Calculate  $log P(a^T \mid s^T, \pi)$
real PolicyRewardBelief::logLikelihood(const Demonstrations<int, int>&D,
						   const  DiscretePolicy& policy) const
{
	real log_prod = 0.0;
	for (uint i = 0; i != D.trajectories.size(); ++i) {
		const Trajectory<int, int>& trajectory = D.trajectories[i];
		for (uint t = 0; t != trajectory.size(); ++t) {
			int s = trajectory.state(t);
			int a = trajectory.action(t);
			real p_s_a = policy.getActionProbability(s, a);
			log_prod += log(p_s_a);
		}
	}
	return log_prod;
}

/// Given a reward matrix, sample a new policy.
FixedDiscretePolicy PolicyRewardBelief::samplePolicy(Matrix& R, real beta)
{
	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			mdp.setFixedReward(s, a, R(s,a)); 
		}
	}
	value_iteration.ComputeStateActionValues(epsilon, (int) ceil(log(epsilon * (1 - gamma)) / log(gamma)));
	return FixedSoftmaxPolicy(value_iteration.Q, beta);
}

/// M-H sampler
void PolicyRewardBelief::MHSampler(Demonstrations<int, int>&D, 
								   int n_iterations)
{
	int n_samples = 0;
	real log_likelihood = LOG_ZERO;
	for (int iter=0; iter<n_iterations; ++iter) {
		Matrix reward = reward_prior.sampleMatrix(); 
		real beta = softmax_prior.generate();
		FixedDiscretePolicy policy = samplePolicy(reward, beta);
		real new_log_likelihood = logLikelihood(D, policy);
		real log_accept_probability = new_log_likelihood - log_likelihood;
		if (log(urandom()) < log_accept_probability) {
			rewards.push_back(reward);
			policies.push_back(policy);
			betas.push_back(beta);
			sample_counts.push_back(1.0);
		} else {
			assert(n_samples); // we must have accepted at least one samplex
			sample_counts[n_samples - 1]++;
		}
	}
}
