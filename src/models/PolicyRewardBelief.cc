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
#include "ExponentialDistribution.h"

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
FixedDiscretePolicy* PolicyRewardBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
	
	return NULL;
}


/// Calculate  $log P(a^T \mid s^T, \pi)$
real PolicyRewardBelief::logLikelihood(const Demonstrations<int, int>&D,
						   const FixedDiscretePolicy& policy) const
{
	real log_prod = 0.0;
	for (uint i = 0; i != D.trajectories.size(); ++i) {
		const Trajectory<int, int>& trajectory = D.trajectories[i];
		for (uint t = 0; t != trajectory.size(); ++t) {
			int s = trajectory.state(t);
			int a = trajectory.action(t);
			real p_s_a = policy.getActionProbability(s, a);
            assert(p_s_a > 0);
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
								   int n_iterations, int n_chains)
{
	int n_samples = 0;
	for (int chain=0; chain<n_chains; ++chain) {
		real log_likelihood = -RAND_MAX;
		for (int iter=0; iter<n_iterations; ++iter) {
			Matrix reward = reward_prior.sampleMatrix(); 
			real beta = softmax_prior.generate();
			FixedDiscretePolicy policy = samplePolicy(reward, beta);
			real new_log_likelihood = logLikelihood(D, policy);
			new_log_likelihood += softmax_prior.log_pdf(beta);
			new_log_likelihood += reward_prior.log_pdf(reward);
			real log_accept_probability = new_log_likelihood - log_likelihood;
			real Z = urandom();
			if (log(Z) < log_accept_probability) {
				rewards.push_back(reward);
				policies.push_back(policy);
				betas.push_back(beta);
				sample_counts.push_back(1.0);
				//logmsg ("New likelihood: %f (%f)\n", new_log_likelihood, log_likelihood);
				log_likelihood = new_log_likelihood;
				n_samples++;
			} else {
				assert(n_samples); // we must have accepted at least one samplex
				sample_counts[n_samples - 1]++;
			}
		}
	}
}

/// Monte-Carlo sampler
void PolicyRewardBelief::MonteCarloSampler(Demonstrations<int, int>&D, 
										   int n_iterations)
{
	int n_samples = 0;
	logmsg (" Running Monte Carlo sampler with %d iterations\n", n_iterations);
	for (int iter=0; iter<n_iterations; ++iter) {
		Matrix reward = reward_prior.sampleMatrix(); 
		real beta = softmax_prior.generate();
		FixedDiscretePolicy policy = samplePolicy(reward, beta);
		rewards.push_back(reward);
		policies.push_back(policy);
		betas.push_back(beta);
		real new_log_likelihood = logLikelihood(D, policy);
		sample_counts.push_back(exp(new_log_likelihood));
		n_samples++;
	}
}


FixedDiscretePolicy* PolicyRewardBelief::getPolicy() 
{
    int n_samples = betas.size();
    real sum = 0.0;
    Matrix R(n_states, n_actions);
    assert(n_samples == rewards.size());
    assert(n_samples == sample_counts.size());
    for (int i=0; i<n_samples; ++i) {
        R += rewards[i] * sample_counts[i];
        sum += sample_counts[i];
    }
    R = R / sum;
    
    for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			mdp.setFixedReward(s, a, R(s,a)); 
		}
	}
    value_iteration.ComputeStateActionValues(epsilon, (int) ceil(log(epsilon * (1 - gamma)) / log(gamma)));
    return value_iteration.getPolicy();


}
