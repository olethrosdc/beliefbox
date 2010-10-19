/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "RSAPI.h"
#include "Rollout.h"
#include "MountainCar.h"
#include "RandomPolicy.h"
#include "MersenneTwister.h"

int main(void)
{
	MersenneTwisterRNG rng;

	// Create a new environment
	Environment<Vector, int>* environment;

	environment = new MountainCar();
	
	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;
	
	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), rng);
	
	// Test how rollouts work with a single rollout sate
	RolloutState state(*environment, environment->getState());
	
	// Setup one rollout per action
	for (uint a=0; a<environment->getNActions(); ++a) {
		state.NewRollout(*policy, a);
	}

	// Sample from this state!
	state.Sample(100);


	delete policy;
	delete environment;
}

#endif
