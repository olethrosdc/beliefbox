// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BPSRModel.h"

BPSRModel::BPSRModel  (int n_obs_, int n_actions_, std::vector<real> rewards_, int tree_depth)
    : n_obs(n_obs_),
      n_actions(n_actions_),
      rewards(rewards_)
{
    mdp_dbg("Creating BPSRModel with %d observations, %d actions and %d rewards. Tree depth %d\n",  n_obs, n_actions, rewards.size(), tree_depth);
    n_rewards = rewards.size();
    std::vector<int> sizes(2);
    sizes[0] = n_obs;
    sizes[1] = n_rewards;
    Z = new DiscreteVector(sizes);
	printf("Making new BPSR with %d observations, %d rewards, %d compound observations, %d actions\n",
		   n_obs, n_rewards, Z->getNCombinations(), n_actions);
    bpsr = new BayesianPredictiveStateRepresentation(Z->getNCombinations(), n_actions, tree_depth, 0.5);
}

BPSRModel::~BPSRModel()
{
    delete bpsr;
    delete Z;
}


/** Observe the initial observation and reward.

    That is, observe \f$x_{t+1}, r_{t+1}\f$.

    The result is passed to the BPSR model by discretising the
    $x_{t+1}, r_{t+1}$ vector.
 */
void BPSRModel::Observe(int x, real r)
{
    std::vector<int> z = getIndexVector(x, r);
    bpsr->Observe(Z->getIndex(z));
}

/** Observe action taken at time t, and the resulting observation and
    reward at the next time step.

    That is, observe \f$a_t, x_{t+1}, r_{t+1}\f$.

    The result is passed to the BPSR model by discretising the
    $x_{t+1}, r_{t+1}$ vector.
 */
void BPSRModel::Observe(int a, int x, real r)
{
    std::vector<int> z = getIndexVector(x, r);
    bpsr->Observe(a, Z->getIndex(z));
}


//real BPSRModel::getTransitionProbability(std::vector<int> history, int a, int x, real r) const
//{
//
//}

/** Obtain the probability of a particular observation and reward at
    the next time step, given a specific current action.

    That is, obtain \f$P(x_{t+1}=x, r_{t+1}=r | a_t = a)\f$.

    This is performed by discretising the $x_{t+1}, r_{t+1}$ vector
    and passing it to the BPSR model.
 */
real BPSRModel::getTransitionProbability(int a, int x, real r) const
{
    std::vector<int> z = getIndexVector(x, r);
    real p = bpsr->ObservationProbability(a, Z->getIndex(z));
    printf ("%f\n", p);
    return p;
}

/** Obtain the probability of a particular observation and reward at
    the next time step, given a specific current action.

    That is, obtain \f$P(x_{t+1}=x, r_{t+1}=r | a_t = a)\f$.

    This is performed by discretising the $x_{t+1}, r_{t+1}$ vector
    and passing it to the BPSR model.
 */
real BPSRModel::getExpectedReward (int a) const
{
    real Er = 0;
    for (int x=0; x<n_obs; ++x) {
        for (uint i=0; i<rewards.size(); ++i) {
            real P = getTransitionProbability(a, x, rewards[i]);
            Er += P * rewards[i];
        }
    }
    return Er;
}

void BPSRModel::Reset()
{
    bpsr->Reset();
}


