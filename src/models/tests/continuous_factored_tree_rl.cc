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

/** \file 
    
    This program discretises a continuous MDP with state space \f$S\f$
    via a projection \f$D : S \to X\f$. Then the problem becomes a POMDP
    with state space \f$S\f$ and observation space \f$X\f$.

    It then uses a context-tree-based POMDP model to perform
    predictions and q-value-function approximation.
 */
#ifdef MAKE_MAIN

#include "FactoredPredictorRL.h"
#include "Random.h"
#include "MersenneTwister.h"
#include "ContinuousStateContextTreeRL.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "Environment.h"

#include <cstdlib>
#include <cstdio>
#include <string>

struct Statistics
{

    std::vector<real> probability;
    std::vector<real> error;
    std::vector<real> reward;

    Statistics(int T) : probability(T), error(T), reward(T)
    {
        for (int t=0; t<T; ++t) {
            probability[t] = 0;
            error[t] = 0;
            reward[t] = 0;
        }
    }
};


bool EvaluateGeneral(Environment<Vector, int>& environment,
                     real action_randomness,
                     FactoredPredictorRL<Vector, int>* factored_predictor,
                     Statistics& statistics);

enum EnvironmentType {
    MOUNTAIN_CAR,
    PENDULUM
};

int main(int argc, char** argv)
{

    int n_actions = 4;
    int n_obs;
	
    if (argc <= 7) {
        fprintf (stderr, "Usage: factored_model environment model_name n_iter T depth env_rand act_rand [environment parameters]\n\
  - environment in {MountainCar, Pendulum}\n\
  - model_name in {BVMM}\n\
  - environment parameters:\n");
        return -argc;
    }
	
	for (int i=0; i<argc; ++i) {
		printf("%d: %s\n", i, argv[i]);
	}
    std::string environment_name(argv[1]);
    std::string model_name(argv[2]);

    int n_iter = atoi(argv[3]);
    int T = atoi(argv[4]);
    int max_depth = atoi(argv[5]);
    real maze_random = atof(argv[6]);
    real action_random = atof(argv[7]);

    EnvironmentType environment_type;
    
    if (!environment_name.compare("MountainCar")) {
        environment_type = MOUNTAIN_CAR;
    } else if (!environment_name.compare("Pendulum")) {
        environment_type = PENDULUM;
    } else {
        Serror ("Unknown environment %s\n",environment_name.c_str());
        exit(-1);
    }

    MersenneTwisterRNG mersenne_twister;
    Statistics statistics(T);
    
    //real prior=0.5;

    for (int iter=0; iter<n_iter; ++iter) {
        mersenne_twister.manualSeed(true_random_bits(false));
        MountainCar mountain_car;
        Pendulum pendulum;
		Vector L_SA;
		Vector U_SA;
		Vector L_S;
		Vector U_S;
		Vector L_A;
		Vector U_A;
        if (environment_type == MOUNTAIN_CAR) {
            n_obs = mountain_car.getNStates();
            n_actions = mountain_car.getNActions();
			L_SA = mountain_car.StateActionLowerBound();
			U_SA = mountain_car.StateActionUpperBound();
			L_S = mountain_car.StateLowerBound();
			U_S = mountain_car.StateUpperBound();
			L_A = mountain_car.ActionLowerBound();
			U_A = mountain_car.ActionUpperBound();
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        } else if (environment_type == PENDULUM) {
            n_obs = pendulum.getNStates();
            n_actions = pendulum.getNActions();
			L_SA = pendulum.StateActionLowerBound();
			U_SA = pendulum.StateActionUpperBound();
			L_S = pendulum.StateLowerBound();
			U_S = pendulum.StateUpperBound();
			L_A = pendulum.ActionLowerBound();
			U_A = pendulum.ActionUpperBound();
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        }
		printf("L_SA: "); L_SA.print(stdout);
		printf("U_SA: "); U_SA.print(stdout);
		printf("L_S: "); L_S.print(stdout);
		printf("U_S: "); U_S.print(stdout);
		printf("L_A: "); L_A.print(stdout);
		printf("U_A: "); U_A.print(stdout);

        FactoredPredictorRL<Vector, int>* factored_predictor; 
        if (!model_name.compare("BVMM")) {
            factored_predictor = new ContinuousStateTFactoredPredictorRL<ContinuousStateContextTreeRL>(n_actions, L_S, U_S, max_depth + 1, max_depth + 1);
        } else {
            fprintf(stderr, "Unrecognised model name %s\n", model_name.c_str());
            exit(-1);
        }
        bool success;
        switch(environment_type) {
        case MOUNTAIN_CAR:
            {

                success = EvaluateGeneral(mountain_car,
                                       action_random,
                                       factored_predictor,
                                       statistics);
            }
            break;
        case PENDULUM:
            {

                success = EvaluateGeneral(pendulum,
                                       action_random,
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
        printf ("%f %f %f #STATS\n",
                statistics.probability[t] * inv_iter,
                statistics.reward[t] * inv_iter, 
                statistics.error[t] * inv_iter);
    }
    return 0;
}


/// This function views looks at a standard MDP-like thing.
bool EvaluateGeneral(Environment<Vector, int>& environment, 
                     real action_randomness,
                     FactoredPredictorRL<Vector, int>* factored_predictor,
                     Statistics& statistics)
{
    int n_obs = environment.getNStates();
    int n_actions = environment.getNActions();
    environment.Reset();
    factored_predictor->Observe(environment.getState());
    int last_action = rand()%n_actions;
    int T = statistics.probability.size();
    printf ("# states: %d, actions: %d, T: %d, a_rand: %f\n", n_obs, n_actions, T, action_randomness);

    std::vector<real> obs_probs(n_obs); // observation probabilities
    std::vector<real> Q(n_actions);

    bool running = false;
    for (int t=0; t<T; ++t) {
        int action = last_action;
        if (!running) {
            factored_predictor->Observe(environment.getState());
            environment.Reset();
            factored_predictor->Reset();
            running = false;
        }

        if (urandom() < action_randomness) {
            action = rand()%n_actions;
		} else {
			for (int a = 0; a < n_actions; ++a) {
				Q[a] = factored_predictor->QValue(a);
			}
			action = ArgMax(Q);
		} 

        running = environment.Act(action);

        Vector observation = environment.getState();
        real reward = environment.getReward();
        //printf ("%d %f ", observation, reward);
        real p = factored_predictor->Observe(action, observation, reward);
        real td_error = factored_predictor->QLearning(0.1, 0.95);
        //assert(p==obs_probs[observation]);
        statistics.probability[t] += p;
        statistics.reward[t] += reward;
        statistics.error[t] += td_error;
        //printf ("%d %d %f\n", action, observation, reward);
    }

    return true;
}

#endif
