/* -*- Mode: C++; -*- */
// copyright (c) 2006-2017 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifdef MAKE_MAIN

/// Other things
#include "MersenneTwister.h"


/// Bayesian RL includes
#include "DiscreteMDPCounts.h"

/// The main thing to test
#include "TreeBRL.h"

/// The basic environment
#include "DiscreteChain.h"

int main(void)
{
	int n_states = 5;
	int n_actions = 2;
	int planning_horizon = 10;
	real discounting = 0.9;

	DiscreteEnvironment* environment = new DiscreteChain(n_states);
	n_states = environment->getNStates();
	n_actions = n_environment->getNActions();
	Matrix rewards(n_states, n_actions);
	for (int s=0; s<n_states; ++s) {
		for (int a=0; s<n_actions; ++a) {
			rewards(s,a) = environment->getExpectedReward(s, a);
		}
	}
	DiscreteMDPCounts belief(n_states, n_actions);

	belief -> setFixedRewards(rewards);

    RandomNumberGenerator* rng;

    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;


	TreeBRL tree (n_states, n_actions, discounting, &belief, rng, planning_horizon);

	

	// Set state to 0
	tree.Reset(0);


	return 0;
}

#endif
