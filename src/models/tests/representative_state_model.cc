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

/** Test that the BPSR model predicts the next observations well.
    
    The BPSRModel gives probabilities of next observations.
    It needs a full history of observations and actions to predict
    the next observation.
    
    However, it always has observations up to x_t and actions
    up to a_{t-1}.

    The BayesianPredictiveStateRepresentation shares the same
    problem. The main difficulty is that the context in the BPSR
    is defined via Factored Markov Chains.

    The context in a FMC is the observation-action history from time
    t-D to time t. It is necessary to have the context in the FMC in
    order to find the right node in the context tree.

 */

#ifdef MAKE_MAIN

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include "Demonstrations.h"
#include "RepresentativeStateModel.h"
#include "Grid.h"
#include "RandomMDP.h"
#include "Gridworld.h"
#include "OneDMaze.h"
#include "DiscreteChain.h"
#include "RiverSwim.h"
#include "OptimisticTask.h"
#include "InventoryManagement.h"
#include "DoubleLoop.h"
#include "ValueIteration.h"
#include "MersenneTwister.h"
#include "DiscreteChain.h"
#include "InventoryManagement.h"

real TestRandomMDP(int n_states, int n_samples);
real TestChainMDP(int n_states, int n_samples);

int main(int argc, char** argv)
{
    
    int n_states = 32;
    int n_samples = 16;
    int n_runs = 10;
    
    logmsg("Usage: representative_state_model [n_states [n_samples [n_runs]]]\n");
    if (argc>1) {
        n_states = atoi(argv[1]);
    }
    if (argc>2) {
        n_samples = atoi(argv[2]);
    }
    if (argc>3) {
        n_runs = atoi(argv[3]);
    }
    Vector ChainMDP_error(n_runs);
    Vector RandomMDP_error(n_runs);
    for (int i=0; i<n_runs; ++i) {
        ChainMDP_error(i) = TestChainMDP(n_states, n_samples);
        RandomMDP_error(i) = TestRandomMDP(n_states, n_samples);
    }

    for (int i=0; i<n_runs; ++i) {
        printf ("%f ", ChainMDP_error(i));
    }
    printf("# chain error\n");

    for (int i=0; i<n_runs; ++i) {
        printf ("%f # chain error\n", RandomMDP_error(i));
    }
    printf("# random MDP error\n");
    
    printf ("%f %f # mean error\n",
            ChainMDP_error.Sum() / (real) n_runs,
            RandomMDP_error.Sum() / (real) n_runs);

    printf("\n\n");
    return 0;
}

real TestRandomMDP(int n_states, int n_samples)
{
	uint n_actions = 2;
	real randomness = 0.1;
	real step_value = -0.1;
	real pit_value = -1;
	real goal_value = 1.0;
	real discount_factor = 0.95;
    real accuracy = 1e-9;
	MersenneTwisterRNG rng;
    RandomMDP random_mdp(n_actions,
                         n_states,
                         randomness,
                         step_value,
                         pit_value,
                         goal_value,
                         &rng);
    Gridworld gridworld("/home/olethros/projects/beliefbox/dat/maze1",
                        randomness,
                        pit_value,
                        goal_value,
                        step_value);

    InventoryManagement inventory(n_states, 32, 0.1, 0.1);

    real total_error = 0;

    DiscreteEnvironment& environment = gridworld;
    printf("%f\n", gridworld.getExpectedReward(0, 0));
    printf("%f\n", environment.getExpectedReward(0, 0));
    n_states = environment.getNStates();
    DiscreteMDP* mdp = environment.getMDP();

    logmsg("Creating model");
#if 0
    // use this if you want to use the actual MDP mdoel
    RepresentativeStateModel<DiscreteMDP, int, int> representative_model(discount_factor, accuracy, *mdp, n_samples, n_actions);
#else
    // use this for sampling
    RepresentativeStateModel<DiscreteEnvironment, int, int> representative_model(discount_factor, accuracy, environment, n_samples, n_actions);
#endif

#if 1

    logmsg("Computing values");
    ValueIteration VI(mdp, discount_factor);
    VI.ComputeStateValuesStandard(accuracy);

    representative_model.ComputeStateValues();

    for (int i=0; i<n_states; ++i) {
        real V = VI.getValue(i);
        real V_approx = representative_model.getValue(i);
        printf ("%f %f # state-value\n", V, V_approx);
        total_error += abs(V - V_approx);
    }

    logmsg("state-action values\n");
    for (int i=0; i<n_states; ++i) {
        for (int a=0; a<n_actions; ++a) {
            real V = VI.getValue(i, a);
            real V_approx = representative_model.getValue(i, a);
            printf ("%d %d %f %f %f\n", i, a, V, V_approx, mdp->getExpectedReward(i, a));
        }
    }
#endif

    delete mdp;
    return total_error;
}



real TestChainMDP(int n_states, int n_samples)
{
	real discount_factor = 0.95;
    real accuracy = 1e-6;
	MersenneTwisterRNG rng;

    logmsg("Creating chain\n");
    DiscreteChain chain(n_states);
    DiscreteEnvironment& environment = chain;

    logmsg("Building model\n");
    DiscreteMDP* mdp = chain.getMDP();
    DiscreteMDP& rmdp = *mdp;
    logmsg("Creating Representative States\n");
    RepresentativeStateModel<DiscreteMDP, int, int> representative_model(discount_factor, accuracy, rmdp, (uint) n_samples, (uint) environment.getNActions());
    
    logmsg("Value iteration\n"); fflush(stdout);
    ValueIteration VI(mdp, discount_factor);
    //VI.ComputeStateValues(accuracy);
    VI.ComputeStateValuesStandard(accuracy);

    logmsg("Approximate VI\n");
    representative_model.ComputeStateValues();

    real total_error = 0;
    logmsg("state values\n");
    for (int i=0; i<n_states; ++i) {
        real V = VI.getValue(i);
        real V_approx = representative_model.getValue(i);
        printf ("%f %f # state-value\n", V, V_approx);
        total_error += abs(V - V_approx);
    }

    logmsg("state-action values\n");
    for (int i=0; i<n_states; ++i) {
        for (int a=0; a<(int) environment.getNActions(); ++a) {
            real V = VI.getValue(i, a);
            real V_approx = representative_model.getValue(i, a);
            printf ("%d %d %f %f %f\n", i, a, V, V_approx, mdp->getExpectedReward(i, a));
        }
    }


    delete mdp;
    return total_error;
}

#endif
