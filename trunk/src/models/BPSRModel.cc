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
        std::vector<int> sizes(3);
        sizes[0] = n_obs;
        sizes[1] = n_actions;
        sizes[2] = n_rewards;
        Z = new DiscreteVector(sizes);
        bpsr = new BayesianPredictiveStateRepresentation(Z->size(), tree_depth, 0.5, false);
    }

BPSRModel::~BPSRModel()
{
    delete bpsr;
    delete Z;
}


void BPSRModel::AddTransition(int a, int x, real r)
{
    std::vector<int> z = getIndex(a, x, r);
    bpsr->Observe(Z->getIndex(z));
}

real BPSRModel::getTransitionProbability(std::vector<int> history, int a, int x) const
{
    return 0;
}
real BPSRModel::getExpectedReward (std::vector<int> history) const
{
    return 0;
}
void BPSRModel::Reset()
{
}
