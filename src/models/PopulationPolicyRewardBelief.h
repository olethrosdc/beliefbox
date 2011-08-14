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
#ifndef POPULATION_POLICY_REWARD_BELIEF_H
#define POPULATION_POLICY_REWARD_BELIEF_H

#include "DiscreteMDP.h"
#include "Demonstrations.h"
#include "DiscretePolicy.h"
#include "RewardDistribution.h"
#include "RewardBelief.h"
#include "DiscretePolicy.h"
#include "ValueIteration.h"

#include <vector>

class DiscreteSpaceRewardDistribution;

/** Prior on policy given reward belief.

 */
class PopulationPolicyRewardBelief
{
protected:
    DiscreteMDP mdp; ///< the actual MDP (transitions assumed known here)
    int n_states; ///< the number of states
    int n_actions; ///< the number of actions
    real eta; ///< the hyperprior parameter (only one!)
    real gamma; ///< value of gamma (assumed known here)
    real epsilon; ///< accuracy
    ValueIteration value_iteration;

    // -- the following values are stored during the sampling procedure -- //
    // samples for each iteration
    std::vector<real> lambdas; ///< values of lambda sampled


    // samples for each demonstration //
    std::vector<std::vector<Matrix> > rewards; ///< matrix of reward function samples ordered [task][sample]
    std::vector<std::vector<FixedDiscretePolicy> > policies; ///< storage for sampled policies from the belief
    std::vector<std::vector<real> > betas; ///< values of beta sampled
    Matrix log_P; ///< log-weight of each sample-demonstration pair
    Vector log_q; ///<  log-weight of each sample

public:
    PopulationPolicyRewardBelief(real eta_,
                                 real gamma_,
                                 const DiscreteMDP& mdp_);
    virtual ~PopulationPolicyRewardBelief();
    
    virtual real logLikelihood(const Demonstrations<int, int>&D,
                               const  FixedDiscretePolicy& policy) const;

    /// Calculate  $log P(a^T \mid s^T, \pi)$
    virtual real Likelihood(const Demonstrations<int, int>&D,
                            const  FixedDiscretePolicy& policy) const
    {
        return exp(logLikelihood(D, policy));
    }

    virtual FixedDiscretePolicy samplePolicy(Matrix& R, real beta);
    virtual std::vector<FixedDiscretePolicy*> getPolicy();
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

#endif
