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
                             RandomNumberGenerator* rng_)
	: ContinuousStateEnvironment(n_features_, n_actions_),
	  mean(n_actions),
	  var(n_actions),
      G(n_states, n_states),
      rng(rng_),
      U_x(n_states),
      L_x(n_states)
{ 
    state.Resize(n_states);

    for (uint i=0; i<n_states; ++i) {
        U_x(i) = 100;
        L_x(i) = -100;
    }
    // this bandit is binary
	for (uint j=0; j<n_actions; ++j) {
		mean[j].Resize(n_states);
		for (uint i=0; i<n_states; ++i) {
			mean[j](i) = normal.generate();
		}
	}

	GammaDistribution gamma(1,1);
	for (uint j=0; j<n_actions; ++j) {
		var(j) = gamma.generate();
	}

    for (uint i=0; i<n_states; ++i) {
        for (uint j=i; j<n_states; ++j) {
            real X = rng->uniform();
            G(i,j) = X;
            G(j,i) = X;
        }
    }
    GenerateContext();
}

LinearContextBandit::~LinearContextBandit()
{
}

void LinearContextBandit::GenerateContext()
{
    Vector y(n_states);
    for (uint j=0; j<n_states; ++j) {
        y(j) = urandom(-1,1);
    }
    const Matrix& rG = G;
    const Vector& ry = y;
    state = rG * ry;
}

/// put the environment in its natural state
void LinearContextBandit::Reset()
{
    GenerateContext();
    reward = 0;
}

/// returns true if the action succeeds, false if we are in a terminal state
bool LinearContextBandit::Act(int action)
{
    normal.setVariance(var(action));
    normal.setMean(Product(&mean[action], &state));
    
    reward = normal.generate();
    //printf("reward: %f %f\n", reward, Product(&mean[action], &state));
    GenerateContext();
    
    return true;  // we continue
}
