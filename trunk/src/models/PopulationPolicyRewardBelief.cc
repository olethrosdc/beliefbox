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

#if 1
#include "PopulationPolicyRewardBelief.h"
#include "ExponentialDistribution.h"


/// Create the belief
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
    Serror("Not implemented\n");
    exit(-1);
    return;
}

/// Monte-Carlo sampler
void PopulationPolicyRewardBelief::MonteCarloSampler(Demonstrations<int, int>&D, 
                                                     int n_iterations)
{
    int n_demonstrations = D.size();

    rewards.resize(n_iterations);
    policies.resize(n_iterations);
    betas.resize(n_iterations);

    log_P.Resize(n_iterations, n_demonstrations);
    log_q.Resize(n_iterations);
    
    logmsg (" Running multi-task Monte Carlo sampler with %d iterations\n", n_iterations);
    ExponentialDistribution softmax_hyperprior(eta);
    ExponentialDistribution reward_hyperprior(1.0);
    for (int iter=0; iter<n_iterations; ++iter) {
        real lambda = softmax_hyperprior.generate();
        ExponentialDistribution softmax_prior(lambda);
        Vector reward_prior_parameters (n_states * n_actions);
        for (int i=0; i<reward_prior_parameters.Size(); ++i) {
            reward_prior_parameters(i) = reward_hyperprior.generate();
        }
        DirichletRewardBelief reward_prior(n_states, n_actions, reward_prior_parameters);

        for (int d=0; d<n_demonstrations; ++d) {
            Matrix reward = reward_prior.sampleMatrix(); 
            real beta = softmax_prior.generate();
            FixedDiscretePolicy policy = samplePolicy(reward, beta);
            rewards[iter].push_back(reward);
            policies[iter].push_back(policy);
            betas[iter].push_back(beta);
            real new_log_likelihood = logLikelihood(D, policy);
            log_P(iter, d) = new_log_likelihood; // log (p_m^k)
        }
    }

    for (int iter=0; iter<n_iterations; ++iter) {
        log_q(iter) = 0;
        for (int d=0; d<n_demonstrations; ++d) {
            log_q(iter) += log_P(iter, d);
        }
    }
    log_q -= log_q.logSum();
}


std::vector<FixedDiscretePolicy*> PopulationPolicyRewardBelief::getPolicy() 
{
    int n_iterations = log_P.Rows();
    int n_demonstrations = log_P.Columns();

    std::vector<FixedDiscretePolicy*> policy_vector;

    for (int d=0; d<n_demonstrations; ++d) {
        // store the average reward in this matrix
        Matrix R(n_states, n_actions);
        for (int i=0; i<n_iterations; ++i) {
            R += rewards[i][d] * exp(log_q(i));
        }
        // solve mean MDP
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                mdp.setFixedReward(s, a, R(s,a)); 
            }
        }
        // find optimal policy for mean MDP
        value_iteration.ComputeStateActionValues(epsilon, (int) ceil(log(epsilon * (1 - gamma)) / log(gamma)));
        // save policy
        policy_vector.push_back(value_iteration.getPolicy());
    }


}

#endif
