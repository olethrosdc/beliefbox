/* -*- Mode: C++; -*- */
// copyright (c) 2017-2020 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
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

/// The basic environments
#include "ContextBandit.h"
#include "DiscreteChain.h"
#include "Blackjack.h"
#include "InventoryManagement.h"
#include "RandomMDP.h"

/// STD
#include <iostream>
#include <memory>
using namespace std;

real RunExperiment(shared_ptr<DiscreteEnvironment> environment,
                   TreeBRL& tree, 
                   int n_steps);

int main(int argc, char** argv) {
    // use a high-quality RNG for the main program
    RandomNumberGenerator* rng;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;

    // use some arbitrary sequence (e.g. a fixed file) for generating environments, to ensure consistency across runs
    DefaultRandomNumberGenerator default_rng;
    RandomNumberGenerator* env_rng = (RandomNumberGenerator*) &default_rng;
    rng->seed();
    env_rng->manualSeed(982374523);
    int n_states = 5;
    int n_actions = 2;
    real discounting = 0.99;
    int n_steps = 1000;

	int state_samples = n_states; ///< number of state samples when branching
    int policy_samples = 2; ///< number of policy/MDP samples at leaf nodes
	int reward_samples = 1; ///< number of reward samples

    // ---- user options ---- //
    int planning_horizon = 2; 
    int leaf_value = 3; 
    int n_experiments = 1; 

	if (argc > 1) {
		planning_horizon = atoi(argv[1]);
	}
	if (argc > 2) {
		leaf_value = atoi(argv[2]);
	}
	if (argc > 3) {
		n_experiments = atoi(argv[3]);
	}
	if (argc > 4) {
		n_steps = atoi(argv[4]);
	}
	if (argc > 5) {
		discounting = atof(argv[5]);
	}
	
    printf("# Making environment\n");
    shared_ptr<DiscreteEnvironment> environment;
    environment = make_shared<DiscreteChain>(n_states, 0);
    //environment = make_shared<ContextBandit>(n_states, n_actions, env_rng, false);
    //environment = make_shared<Blackjack>(env_rng);
    //environment = make_shared<RandomMDP>(n_states, n_actions, 0.1, -0.1, -1, 1, env_rng);
	// environment = make_shared<InventoryManagement>(10, 5, 0.2, 0.1);-=  n_states = environment->getNStates();
    n_actions = environment->getNActions();
	n_states = environment->getNStates();
#if 0
    //simplify things by fixing the rewards
    printf("# Setting up belief\n");
    Matrix rewards(n_states, n_actions);
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            rewards(s,a) = environment->getExpectedReward(s, a);
        }
    }
    belief.setFixedRewards(rewards);
#endif

    printf("# full sampling\n");


	
    real dirichlet_mass = 0.5;
    enum DiscreteMDPCounts::RewardFamily reward_prior = DiscreteMDPCounts::BETA;
    DiscreteMDPCounts belief(n_states, n_actions, dirichlet_mass, reward_prior);
	//NullMDPModel belief(n_states, n_actions);
	
    Vector U(n_experiments);
    for (int experiment=0;
         experiment<n_experiments;
         experiment++) {
        
        TreeBRL tree (n_states, n_actions, discounting, &belief, rng, planning_horizon, state_samples, policy_samples, reward_samples, (TreeBRL::LeafNodeValue) leaf_value);
        // Set state to 0

        real total_reward = RunExperiment(environment, tree, n_steps);
        printf("H:%d,\tV:%d,\tR:%f\n", planning_horizon, leaf_value,total_reward);
        U(experiment) = total_reward;
    }
    printf("L:%f,\tM:%f,\tU:%f\n",
           Min(U),
           U.Sum() / (real) U.Size(),
           Max(U));
    return 0;
}

real RunExperiment(shared_ptr<DiscreteEnvironment> environment,
                   TreeBRL& tree, 
                   int n_steps)
{
    environment->Reset();
    tree.Reset(environment->getState());
    real total_reward = 0;
    for (int t=0; t<n_steps; ++t) {
        int state = environment->getState();
        real reward = environment->getReward();
        int action = tree.Act(reward, state);
        //			action = 1;
        bool action_OK = environment->Act(action);
        total_reward += reward;
        printf("%d %d %f  # reward\n", state, action, reward);
        if (!action_OK) {
            state = environment->getState();
            reward = environment->getReward();
            tree.Act(reward, state);
            environment->Reset();
            tree.Reset(environment->getState());
        }
    }
    return total_reward;
}

#endif
