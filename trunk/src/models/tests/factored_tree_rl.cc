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

#include "FactoredPredictorRL.h"
#include "POMDPGridworld.h"
#include "Corridor.h"
#include "Random.h"
#include "MersenneTwister.h"
#include "DiscretePOMDP.h"
#include "POMDPBeliefState.h"
#include "POMDPBeliefPredictor.h"
#include "ContextTreeRL.h"

#include <cstdlib>
#include <cstdio>
#include <string>

struct Statistics
{

    std::vector<real> probability;
    std::vector<real> error;

    Statistics(int T) : probability(T), error(T)
    {
        for (int t=0; t<T; ++t) {
            probability[t] = 0;
            error[t] = 0;
        }
    }
};

bool EvaluateMaze(std::string maze, 
                  int n_obs,
                  real maze_random,
                  real action_random,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictorRL* factored_predictor,
                  Statistics& statistics);

bool Evaluate1DWorld(Corridor& environment,
                     real action_randomness,
                     FactoredPredictorRL* factored_predictor,
                     Statistics& statistics);

enum EnvironmentType {
    ONEDMAZE,
    GRIDWORLD
};

int main(int argc, char** argv)
{

    int n_actions = 4;
    int n_obs;
	
    if (argc <= 8) {
        fprintf (stderr, "Usage: factored_model environment model_name n_iter T depth env_rand act_rand [environment parameters]\n\
  - environment in {1DMaze, Gridworld}\n\
  - model_name in {FMC, BFMC, CTW, BVMM}\n\
  - environment parameters:\n\
    - 1DMaze: number_of_states\n\
    - Gridworld: maze, n_obs in {2, 16}\n");
        return -argc;
    }
	

    std::string environment_name(argv[1]);
    std::string model_name(argv[2]);

    int n_iter = atoi(argv[3]);
    int T = atoi(argv[4]);
    int max_depth = atoi(argv[5]);
    real maze_random = atof(argv[6]);
    real action_random = atof(argv[7]);

    EnvironmentType environment_type;

    if (!environment_name.compare("1DMaze")) {
        if (argc!=9) {
            Serror("Model parameters: number of states\n");
            exit(-1);
        }
        environment_type = ONEDMAZE;
        n_actions = 2;
        n_obs = 2;
    } else if (!environment_name.compare("Gridworld")) {
        if (argc!=10) {
            Serror("Model parameters: maze_filename, n_obs in {2, 16}\n");
            exit(-1);
        }
        environment_type = GRIDWORLD;
        n_actions = 4;
        n_obs = atoi(argv[9]);
    } else {
        Serror ("Unknown environment %s\n",environment_name.c_str());
        exit(-1);
    }

    MersenneTwisterRNG mersenne_twister;
    Statistics statistics(T);
    
    //real prior=0.5;

    for (int iter=0; iter<n_iter; ++iter) {
        mersenne_twister.manualSeed(true_random_bits(false));

        FactoredPredictorRL* factored_predictor; 
        if (!model_name.compare("BVMM")) {
            factored_predictor = new TFactoredPredictorRL<ContextTreeRL>(n_actions, n_obs, max_depth + 1);
        } else {
            fprintf(stderr, "Unrecognised model name %s\n", model_name.c_str());
            exit(-1);
        }
        bool success;
        switch(environment_type) {
        case ONEDMAZE:
            {
                int n_states = atoi(argv[8]);
                Corridor environment(n_states, maze_random, &mersenne_twister);
                //if (!factored_predictor) {
                //                    factored_predictor = new POMDPBeliefPredictor(environment.getPOMDP());
                //}
                success = Evaluate1DWorld(environment,
                                          action_random,
                                          factored_predictor,
                                          statistics);
            }
            break;
        case GRIDWORLD:
            {
                std::string maze(argv[8]);
                int n_obs = atoi(argv[9]);
                success = EvaluateMaze(maze,
                                       n_obs,
                                       maze_random,
                                       action_random,
                                       &mersenne_twister,
                                       factored_predictor,
                                       statistics);
            }
            break;
        default:
            Serror("Unknown world\n");
            exit(-1);
        }
        if (!success) {
            fprintf(stderr, "Failure in iteration %d\n", iter);
        }
        delete factored_predictor;
    }
	
    real inv_iter = 1.0 / (real) n_iter;
    for (int t=0; t<T; ++t) {
        //printf ("%f %f\n", statistics.probability[t] * inv_iter,
        //statistics.error[t] * inv_iter);
        printf ("%f\n", statistics.probability[t] * inv_iter);
    }
    return 0;
}

bool EvaluateMaze(std::string maze, 
                  int n_obs,
                  real world_randomness,
                  real action_randomness,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictorRL* factored_predictor,
                  Statistics& statistics)
{
    POMDPGridworld environment(environment_rng, maze.c_str(),
                               n_obs, 
                               world_randomness);
    int n_actions = environment.getNActions();
    environment.Reset();
    factored_predictor->Observe(environment.getObservation());
    int last_action = 0;
    int T = statistics.probability.size();
    std::vector<real> obs_probs(n_obs); ///< observation probabilities

	for (int t=0; t<T; ++t) {
            int action = last_action;
            
            // Add it here.

            environment.Act(action);
            
#if 0
            for (int i=0; i<n_obs; ++i) {
                obs_probs[i] = factored_predictor->ObservationProbability(action, i);
            }
#endif
            int observation = environment.getObservation();
            real reward = environment.getReward();
            if ((n_obs == 2 && observation) || urandom() < action_randomness) {
                last_action = rand()%n_actions;
            }
            
            real p = factored_predictor->Observe(action, observation, reward);
            assert(p==obs_probs[observation]);
            statistics.probability[t] += p;
#if 0
            if (ArgMax(obs_probs) != observation) {
                statistics.error[t]++;
            }
#endif
	}

	return true;
}


bool Evaluate1DWorld(Corridor& environment,
                     real action_randomness,
                     FactoredPredictorRL* factored_predictor,
                     Statistics& statistics)
{

    int n_states = environment.getNStates();
    int n_actions = environment.getNActions();
    int n_obs = environment.getNObservations();
    environment.Reset();
    factored_predictor->Observe(environment.getObservation());
    int last_action = rand()%n_actions;
    int T = statistics.probability.size();
    printf ("# states: %d, actions: %d, obs: %d, T: %d, a_rand: %f\n", n_states, n_actions, n_obs, T, action_randomness);

    std::vector<real> obs_probs(n_obs); // observation probabilities
    
    for (int t=0; t<T; ++t) {
        int action = last_action;
        environment.Act(action);
        //for (int i=0; i<n_obs; ++i) {
        //obs_probs[i] = factored_predictor->ObservationProbability(action, i, 0);
            //            printf("# %d %f\n", i, obs_probs[i]);
        //}
        int observation = environment.getObservation();
        real reward = environment.getReward();
        if (urandom() < action_randomness) {
            last_action = rand()%n_actions;
        }

        real p = factored_predictor->Observe(action, observation, reward);
        assert(p==obs_probs[observation]);
        statistics.probability[t] += p;
        if (ArgMax(obs_probs) != observation) {
            statistics.error[t]++;
        }
    }

    return true;
}
#endif
