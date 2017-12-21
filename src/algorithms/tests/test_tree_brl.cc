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
#include "Blackjack.h"
#include "InventoryManagement.h"

/// STD
#include <iostream>
#include <memory>
using namespace std;

int main(void) {
    RandomNumberGenerator* rng;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;
    rng->seed();
    int n_states = 5;
    int n_actions = 2;
    int max_planning_horizon = 10;
    real discounting = 0.95;
	//    int n_samples = 2; ///< number of state samples when branching
    //int n_mdp_samples = 2; ///< number of MDP samples at leaf nodes

    printf("# Making environment\n");
    unique_ptr<DiscreteEnvironment> environment;
    environment = make_unique<DiscreteChain>(n_states);
    //environment = make_unique<Blackjack>(rng);
    //environment = make_unique<InventoryManagement>(10, 5, 0.2, 0.1);
    n_states = environment->getNStates();
    n_actions = environment->getNActions();

    printf("# Setting up belief\n");
    Matrix rewards(n_states, n_actions);
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            rewards(s,a) = environment->getExpectedReward(s, a);
        }
    }

    // simplify things by fixing the rewards
    //belief.setFixedRewards(rewards);


    printf("# full sampling\n");
    for (int planning_horizon=0;
         planning_horizon<max_planning_horizon;
         planning_horizon++) {

		real dirichlet_mass = 0.5;
		enum DiscreteMDPCounts::RewardFamily reward_prior = DiscreteMDPCounts::BETA;
		DiscreteMDPCounts belief(n_states, n_actions, dirichlet_mass, reward_prior);
			
		
        TreeBRL tree (n_states, n_actions, discounting, &belief, rng, planning_horizon);
        // Set state to 0
        tree.Reset(0);

		real total_reward = 0;
		for (int t=0; t<10000; ++t) {
			int state = environment->getState();
			real reward = environment->getReward();
			int action = tree.Act(reward, state);
			//			action = 1;
			environment->Act(action);
			total_reward += reward;
			printf("%d %d %f  # reward\n", state, action, reward);
		}
		printf("H: %d, R: %f\n", planning_horizon, total_reward);
    }

    return 0;
}

#endif
