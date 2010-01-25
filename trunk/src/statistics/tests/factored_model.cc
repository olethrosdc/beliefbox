/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
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

#include "FactoredPredictor.h"
#include "FactoredMarkovChain.h"
#include "BayesianFMC.h"
#include "BayesianPredictiveStateRepresentation.h"
#include "BayesianPredictiveStateRepresentationCTW.h"
#include "POMDPGridworld.h"
#include "Random.h"
#include "MersenneTwister.h"

#include <cstdlib>
#include <cstdio>
#include <string>

class Corridor
{
protected:
    int n_states;
    int observation;  ///< current observation
    int state; ///< current state
public:
    Corridor(int n_states_) : n_states(n_states_)
    {
        assert(n_states > 0);
    }
    void Reset()
    {
        state = 0;
        observation = 0;
    }
    int getObservation()
    {
        return observation;
    }
    bool Act(int action)
    {
        switch(action) {
        case 0:
            state--; 
            break;
        case 1:
            state++;
            break;
        default:
            break;
        }
        observation = 0;
        if (state < 0) {
            state = 0;
            observation = 1;
        } else if (state >= n_states) {
            state = n_states - 1;
            observation = 1;
        }
        return 0;
    }
};

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
                  int maze_height,
                  int maze_width,
                  int n_obs,
                  real random,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictor* factored_predictor,
                  Statistics& statistics);

int main(int argc, char** argv)
{

    int n_actions = 4;
    int n_obs = 2;
	
    if (argc != 6) {
        fprintf (stderr, "Usage: factored_model n_iter T states depth model_name\n  - model_name in {FMC, BFMC, CTW, BVMM}\n");
        return -argc;
    }
	
    int n_iter = atoi(argv[1]);
    int T = atoi(argv[2]);
    int n_states = atoi(argv[3]);
    int max_depth = atoi(argv[4]);
    std::string model_name(argv[5]);

    //Corridor environment(n_states);
    //int n_actions = 2;
#if 1
    std::string homedir(getenv("HOME"));
    std::string maze = homedir + "/projects/beliefbox/dat/maze1c";
    int arg_maze_height = 8;
    int arg_maze_width = 8;
    real random = 0.01;//atof(argv[2]);
#else
    std::string maze = argv[4];
#endif
    MersenneTwisterRNG mersenne_twister;
    Statistics statistics(T);
    
    real prior=0.1;

    for (int iter=0; iter<n_iter; ++iter) {
        FactoredPredictor* factored_predictor; 
        if (!model_name.compare("FMC")) {
            factored_predictor = new FactoredMarkovChain(n_actions, n_obs, max_depth);
        } else if (!model_name.compare("BFMC")) {
            factored_predictor = new BayesianFMC(n_obs, n_actions, max_depth + 1, prior);
        } else if (!model_name.compare("BVMM")) {
            factored_predictor = new BayesianPredictiveStateRepresentation(n_obs, n_actions,  max_depth + 1, prior);
        } else if (!model_name.compare("CTW")) {
            factored_predictor = new BayesianPredictiveStateRepresentationCTW(n_obs, n_actions,  max_depth + 1, prior);
        } else {
            fprintf(stderr, "Unrecognised model name %s\n", model_name.c_str());
            exit(-1);
        }
        bool success = EvaluateMaze(maze, arg_maze_height, arg_maze_width,
                                    n_obs,	random, &mersenne_twister,
                                    factored_predictor, statistics);
        if (!success) {
            fprintf(stderr, "Failure in iteration %d\n", iter);
        }
        delete factored_predictor;
    }
	
    real inv_iter = 1.0 / (real) n_iter;
    for (int t=0; t<T; ++t) {
        printf ("%f %f\n", statistics.probability[t] * inv_iter,
                statistics.error[t] * inv_iter);
    }
    return 0;
}

bool EvaluateMaze(std::string maze, 
                  int maze_height,
                  int maze_width,
                  int n_obs,
                  real random,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictor* factored_predictor,
                  Statistics& statistics)
{
    POMDPGridworld environment(environment_rng, maze.c_str(),
                               n_obs, maze_height, maze_width, random);
    int n_actions = environment.getNActions();
    environment.Reset();
    factored_predictor->Observe(environment.getObservation());
    int last_action = 0;
    int T = statistics.probability.size();
    std::vector<real> obs_probs(n_obs); ///< observation probabilities
	for (int t=0; t<T; ++t) {
            int action = last_action;
            environment.Act(action);
            for (int i=0; i<n_obs; ++i) {
                obs_probs[i] = factored_predictor->ObservationProbability(action, i);
            }
            int observation = environment.getObservation();
            if (observation || urandom() < 0.1) {
                last_action = rand()%n_actions;
            }
            real p = factored_predictor->Observe(action, observation);
            assert(p==obs_probs[observation]);
            statistics.probability[t] += p;
            if (ArgMax(obs_probs) != observation) {
                statistics.error[t]++;
            }
	}

	return true;
}
#endif
