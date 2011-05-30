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

#include "RewardBelief.h"
#include "RewardDistribution.h"
#include "Matrix.h"

DiscreteSpaceRewardDistribution* DirichletRewardBelief::sample() const
{
	DiscreteSpaceRewardDistribution* reward_distribution
		= new DiscreteSpaceRewardDistribution (n_states, n_actions);
	Matrix M = sampleMatrix();
	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			reward_distribution->setFixedReward(s, a, M(s, a));
		}
	}
	return reward_distribution;
}

Matrix DirichletRewardBelief::sampleMatrix() const
{

	Matrix M(n_states, n_actions);
	Vector R = dirichlet.generate();
	int i=0; 
	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			M(s, a) = R(i++);
		}
	}
	return M;
}
