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
#include "BPSRModel.h"
#include "POMDPGridworld.h"
#include "MersenneTwister.h"

int main(int argc, char** argv)
{
    std::vector<real> rewards(4);
    rewards[0] = -1.0;
    rewards[1] = -0.1;
    rewards[2] = 0.0;
    rewards[3] = 1.0;
    int n_actions = 4;
    real random = 0.01;
    int tree_depth = atoi(argv[1]);
    std::string homedir(getenv("HOME"));
    std::string maze = homedir + "/projects/beliefbox/dat/maze1";
    POMDPGridworld environment(maze.c_str(), 8, 8, n_actions, random);
    
    int n_obs = environment.getNObs();
    BPSRModel model(n_obs, n_actions, rewards, tree_depth);
    MersenneTwisterRNG mersenne_twister;
    RandomNumberGenerator& rng = mersenne_twister;


    rng.seed();
    int T = 100000;
    int observation = environment.getObservation();
    model.Observe(observation, 0.0);
    for (int t=0; t<T; ++t) {
        //        environment.Show();
        int state = environment.getState();
        int action = rng.discrete_uniform(n_actions);
        bool terminate = environment.Act(action);
        observation = environment.getObservation();
        real reward = environment.getReward();

        real probability = model.getTransitionProbability(action, observation, reward);
        printf ("%d %d %d %f %f # LOG\n", state, observation, action, reward, probability);
        model.Observe(action, observation, reward);
        
        if (terminate) {
            environment.Reset();
        }
     }

}


#endif
