// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "LinearContextBandit.h"
#include "GammaDistribution.h"

LinearContextBandit::LinearContextBandit(uint n_actions_,
                             uint n_features_,
                             RandomNumberGenerator* rng)
	: n_actions(n_actions_),
	  n_features(n_features_),
	  mean(n_actions),
	  std(n_actions)
{ 
    this->rng = rng;
    // this bandit is binary
	NormalDistribution normal(0,1);
	for (uint j=0; j<n_actions; ++j) {
		mean[j].Resize(n_features);
		for (uint i=0; i<n_features; ++i) {
			mean[i](j) = normal.generate();
		}
	}

	GammaDistribution gamma(1,1);
	for (uint j=0; j<n_actions; ++j) {
		std(j) = gamma.generate();
	}
	
}


/// put the environment in its natural state
void LinearContextBandit::Reset()
{
    state = (int) rng->discrete_uniform(n_states);
    reward = 0;
}

/// returns true if the action succeeds, false if we are in a terminal state
bool LinearContextBandit::Act(int action)
{
    real sigma = 2.0;
    normal.setVariance();
    normal.setMean(getMean(action));

    reward = normal.generate();
    //printf("reward: %f\n", reward);
    state = (int) rng->discrete_uniform(n_states);


    return true;  // we continue
}
