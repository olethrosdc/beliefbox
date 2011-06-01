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

#include "PopulationPolicyRewardBelief.h"



/// Create teh blief
PopulationPolicyRewardBelief::PopulationPolicyRewardBelief(real eta_,
									   real gamma_,
									   const DiscreteMDP& mdp_)
	: mdp(mdp_),
	  n_states(mdp.getNStates()),
	  n_actions(mdp.getNActions()),
      eta(eta_),
	  gamma(gamma_),
	  value_iteration(&mdp, gamma)
{
	setAccuracy(1e-3);
}

/// Destroy the belief
PopulationPolicyRewardBelief::~PopulationPolicyRewardBelief()
{
}

/// This samples the posterior belief given the original sample....
FixedDiscretePolicy* PopulationPolicyRewardBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
	
	return NULL;
}


/// Calculate  $log P(a^T \mid s^T, \pi)$
real PopulationPolicyRewardBelief::logLikelihood(const Demonstrations<int, int>&D,
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
FixedDiscretePolicy PopulationPolicyRewardBelief::samplePolicy(Matrix& R, real beta)
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
void PopulationPolicyRewardBelief::MHSampler(Demonstrations<int, int>&D, 
								   int n_iterations, int n_chains)
{
    return;
}

/// Monte-Carlo sampler
void PopulationPolicyRewardBelief::MonteCarloSampler(Demonstrations<int, int>&D, 
										   int n_iterations)
{
	int n_samples = 0;
	logmsg (" Running Monte Carlo sampler with %d iterations\n", n_iterations);
    DirichletRewardBelief reward_hyperprior(n_states, n_actions);
    ExponentialDistribution softmax_hyperprior(eta);
	for (int iter=0; iter<n_iterations; ++iter) {
        real lambda = softmax_hyperprior.generate();
        ExponentialDistribution softmax_prior(lambda);
        DirichletRewardBelief reward_prior(n_states, n_actions,
                                           reward_hyperprior.generate());
        int n_demonstrations = D.trajectories.size();
        for (int d=0; d<n_demonstrations; ++d) {
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
}


FixedDiscretePolicy* PopulationPolicyRewardBelief::getPolicy() 
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
