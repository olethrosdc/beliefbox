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
#include "DiscretePOMDP.h"
#include "POMDPBeliefState.h"
#include "POMDPBeliefPredictor.h"

#include <cstdlib>
#include <cstdio>
#include <string>

class Corridor
{
protected:
    int n_states;
    int observation;  ///< current observation
    int state; ///< current state
    real randomness;
    RandomNumberGenerator* rng;
    DiscretePOMDP* pomdp;
public:
    Corridor(int n_states_, real randomness_, RandomNumberGenerator* rng_)
        : n_states(n_states_),
          randomness(randomness_),
          rng(rng_)
    {
        printf ("# Making Corridor of length %d\n", n_states);
        assert(n_states > 0);
        pomdp = new DiscretePOMDP(n_states, 2, 2);
        int n_obs = 2;
        // clear all
        for (int s=0; s<n_states; ++s) {
            for (int s2=0; s2<n_states; ++s2)  {
                pomdp->setNextStateProbability(s, 0, s2, 0);
                pomdp->setNextStateProbability(s, 1, s2, 0);
            }
            for (int x=0; x<n_obs; ++x)  {
                pomdp->setObservationProbability(s, 0, x, 0);
                pomdp->setObservationProbability(s, 0, x, 0);
            }
        }
        
        // set states
        for (int s=0; s<n_states-1; ++s) {
            pomdp->setNextStateProbability(s+1, 0, s, 1);
            pomdp->setNextStateProbability(s, 1, s+1, 1);
        }
        pomdp->setNextStateProbability(0, 0, 0, 1);
        pomdp->setNextStateProbability(n_states-1, 1, n_states-1, 1);

        // set obs
        for (int s=0; s<n_states-1; ++s) {
            // see nothing when no obstacle
            pomdp->setObservationProbability(s+1, 0, 0, 1 - randomness);
            pomdp->setObservationProbability(s,   1, 0, 1 - randomness);
            pomdp->setObservationProbability(s+1, 0, 1, randomness);
            pomdp->setObservationProbability(s,   1, 1, randomness);
        }
        // see something when obstacle
        pomdp->setObservationProbability(0, 0, 1, 1 - randomness);
        pomdp->setObservationProbability(0, 0, 0, randomness);
        pomdp->setObservationProbability(n_states-1, 1, 1, 1 - randomness);
        pomdp->setObservationProbability(n_states-1,  1, 0, randomness);

        pomdp->check();
    }
    ~Corridor()
    {
        delete pomdp;
    }
    void Reset()
    {
        state = 0;
        observation = 0;
    }
    DiscretePOMDP* getPOMDP()
    {
        return pomdp;
    }
    int getObservation()
    {
        return observation;
    }
    int getNStates()
    {
        return n_states;
    }
    int getNActions()
    {
        return 2;
    }
    int getNObservations()
    {
        return 2;
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
        if (rng->uniform() < randomness) {
            observation = 1 - observation;
        }
        pomdp->setObservation(observation);
        pomdp->setReward(0);
        pomdp->setState(state);
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
                  int n_obs,
                  real maze_random,
                  real action_random,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictor* factored_predictor,
                  Statistics& statistics);

bool Evaluate1DWorld(Corridor& environment,
                     real action_randomness,
                     FactoredPredictor* factored_predictor,
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
    
    real prior=0.5;

    for (int iter=0; iter<n_iter; ++iter) {
        mersenne_twister.manualSeed(true_random_bits(false));

        FactoredPredictor* factored_predictor; 
        if (!model_name.compare("FMC")) {
            factored_predictor = new FactoredMarkovChain(n_actions, n_obs, max_depth);
        } else if (!model_name.compare("BFMC")) {
            factored_predictor = new BayesianFMC(n_obs, n_actions, max_depth + 1, prior);
        } else if (!model_name.compare("BVMM")) {
            factored_predictor = new BayesianPredictiveStateRepresentation(n_obs, n_actions,  max_depth + 1, prior);
        } else if (!model_name.compare("CTW")) {
            factored_predictor = new BayesianPredictiveStateRepresentationCTW(n_obs, n_actions,  max_depth + 1, prior);
        } else if (!model_name.compare("POMDP")) {
            factored_predictor = NULL;
            // to be made later
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
                if (!factored_predictor) {
                    factored_predictor = new POMDPBeliefPredictor(environment.getPOMDP());
                }
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
        printf ("%f %f\n", statistics.probability[t] * inv_iter,
                statistics.error[t] * inv_iter);
    }
    return 0;
}

bool EvaluateMaze(std::string maze, 
                  int n_obs,
                  real world_randomness,
                  real action_randomness,
                  RandomNumberGenerator* environment_rng,
                  FactoredPredictor* factored_predictor,
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
            environment.Act(action);
            for (int i=0; i<n_obs; ++i) {
                obs_probs[i] = factored_predictor->ObservationProbability(action, i);
            }
            int observation = environment.getObservation();
            if (observation || urandom() < action_randomness) {
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


bool Evaluate1DWorld(Corridor& environment,
                     real action_randomness,
                     FactoredPredictor* factored_predictor,
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
        for (int i=0; i<n_obs; ++i) {
            obs_probs[i] = factored_predictor->ObservationProbability(action, i);
            //            printf("# %d %f\n", i, obs_probs[i]);
        }
        int observation = environment.getObservation();
        if (urandom() < action_randomness) {
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
