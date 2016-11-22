// -*- Mode: c++ -*-
// copyright (c) 2016 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMMDP.h"
#include "Random.h"
#include "RandomNumberGenerator.h"
#include "Dirichlet.h"

#include <cstdio>

int main(void)
{
	// Preliminary setup
	printf("Starting up\n");
	DefaultRandomNumberGenerator rng;
	
	// Create the basic MMDP
	int n_states = 2;
	int n_players = 2;
	int n_actions = 2;
	DiscreteMMDP MMDP(n_players, n_states, n_actions);
	int n_joint_actions = MMDP.getNJointActions();

	printf("Creating MMDP with N:%d, S:%d, A:%d, A^N:%d\n",
		   n_players,
		   n_states,
		   n_actions,
		   n_joint_actions);
	   
	// Create the rewards uniformly in [0,1]
	for (int i=0; i<n_states; ++i) {
		MMDP.setReward(i, rng.uniform());
	}

	// Now, for every state and joint-action pair, generate a random
	// transition probability to the next state. This code is
	// specialised for two player games. We'd need to do a combinatorial thing.
	DirichletDistribution dirichlet(n_states); /// Use a Dirichlet to generate multinomials
	for (int i=0; i<n_states; ++i) {
		DiscreteJointAction action(n_players, n_actions);
		for (int a=0; a<n_actions; ++a) {
			action.set(1, a);
			for (int b=0; b<n_actions; ++b) {
				action.set(1, b);
				Vector P = dirichlet.generate(); // put the multinomial in a Vector class
				for (int j=0; j<n_states; ++j) {
					MMDP.setTransitionProbability(i, action, j, P(j));
				}
			}
		}
	}
	printf("Created transition probabilities\n");

	// 1. Test policy evaluation and backwards induction for when the players have the same model, i.e. optimise for the joint policy.
	
	// 2. Consider the case where each player has a different MMDP belief. 


	
	return 0;
}

