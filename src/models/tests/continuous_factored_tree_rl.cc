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
#include "MountainCar3D.h"
#include "Pendulum.h"
#include "LinearContextBandit.h"
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

struct EpisodeStatistics
{
    std::vector<real> total_reward;
    std::vector<real> utility;
    std::vector<real> steps;
    int n_episodes;
    void AddEpisode(real total_reward_, real utility_, real steps_)
    {
        n_episodes++;
        total_reward.push_back(total_reward_);
        utility.push_back(utility_);
        steps.push_back(steps_);
    }

};


bool EvaluateGeneral(Environment<Vector, int>& environment,
                     real action_randomness,
                     FactoredPredictorRL<Vector, int>* factored_predictor, 
                     real discount_factor,
                     EpisodeStatistics& episode_statistics,
                     Statistics& statistics);

enum EnvironmentType {
    MOUNTAIN_CAR,
    MOUNTAIN_CAR_3D,
    PENDULUM,
    LINEAR_BANDIT
};

int main(int argc, char** argv)
{

    int n_actions = 4;
    int n_obs;
	
    if (argc <= 9) {
        fprintf (stderr, "Usage: factored_model environment model_name n_iter T depth env_rand act_rand weight_factor depth_factor[environment parameters]\n\
  - environment in {MountainCar, Pendulum, LinearBandit, MountainCar3D}\n\
  - env_rand environment randomness in [0,1]\n\
  - weight_factor in (0,1] -- default 0.5\n\
  - depth_factor in [0,1] -- default 1\n\
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
    real environment_randomness = atof(argv[6]);
    real action_random = atof(argv[7]);
    real weight_factor = atof(argv[8]);
    real depth_factor = atof(argv[9]);
    real discount_factor = 0.95;

    EnvironmentType environment_type;
    
    if (!environment_name.compare("MountainCar")) {
        environment_type = MOUNTAIN_CAR;
    } else if (!environment_name.compare("MountainCar3D")) {
        environment_type = MOUNTAIN_CAR_3D;
    } else if (!environment_name.compare("Pendulum")) {
        environment_type = PENDULUM;
    } else if (!environment_name.compare("LinearBandit")) {
        environment_type = LINEAR_BANDIT;
    } else {
        Serror ("Unknown environment %s\n",environment_name.c_str());
        exit(-1);
    }

    MersenneTwisterRNG mersenne_twister;
    Statistics statistics(T);
    std::vector<EpisodeStatistics> episode_statistics(n_iter);
    //real prior=0.5;

    for (int iter=0; iter<n_iter; ++iter) {
        mersenne_twister.manualSeed(true_random_bits(false));
        MountainCar mountain_car;
        MountainCar3D mountain_car_3d;
        Pendulum pendulum;
        LinearContextBandit bandit(4,4, &mersenne_twister);        
        mountain_car.setRandomness(environment_randomness);
        mountain_car_3d.setRandomness(environment_randomness);
        pendulum.setRandomness(environment_randomness);
		Vector L_S;
		Vector U_S;
        if (environment_type == MOUNTAIN_CAR) {
            n_obs = mountain_car.getNStates();
            n_actions = mountain_car.getNActions();
			L_S = mountain_car.StateLowerBound() - 1;
			U_S = mountain_car.StateUpperBound() + 1;
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        } else if (environment_type == MOUNTAIN_CAR_3D) {
            n_obs = mountain_car_3d.getNStates();
            n_actions = mountain_car_3d.getNActions();
			L_S = mountain_car_3d.StateLowerBound() - 1;
			U_S = mountain_car_3d.StateUpperBound() + 1;
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        } else if (environment_type == PENDULUM) {
            n_obs = pendulum.getNStates();
            n_actions = pendulum.getNActions();
			L_S = pendulum.StateLowerBound() - 1;
			U_S = pendulum.StateUpperBound() + 1;
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        } else if (environment_type == LINEAR_BANDIT) {
            n_obs = bandit.getNStates();
            n_actions = bandit.getNActions();
			L_S = bandit.StateLowerBound() - 1;
			U_S = bandit.StateUpperBound() + 1;
            printf ("# Total observations: %d, actions: %d\n", n_obs, n_actions);
        }
		printf("L_S: "); L_S.print(stdout);
		printf("U_S: "); U_S.print(stdout);
        
        FactoredPredictorRL<Vector, int>* factored_predictor; 
        if (!model_name.compare("BVMM")) {
            factored_predictor = new ContinuousStateTFactoredPredictorRL<ContinuousStateContextTreeRL>(n_actions, L_S, U_S, max_depth , 8, depth_factor, weight_factor);
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
                                          discount_factor,
                                          episode_statistics[iter],
                                          statistics);
            }
            break;
        case MOUNTAIN_CAR_3D:
            {

                success = EvaluateGeneral(mountain_car_3d,
                                          action_random,
                                          factored_predictor,
                                          discount_factor,
                                          episode_statistics[iter],
                                          statistics);
            }
            break;
        case PENDULUM:
            {

                success = EvaluateGeneral(pendulum,
                                          action_random,
                                          factored_predictor,
                                          discount_factor,
                                          episode_statistics[iter],
                                          statistics);
            }
            break;
        case LINEAR_BANDIT:
            {

                success = EvaluateGeneral(bandit,
                                          action_random,
                                          factored_predictor,
                                          discount_factor,
                                          episode_statistics[iter],
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
        factored_predictor->Show();
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

    int max_episodes = 0;
    for (int iter=0; iter<n_iter; ++iter) {
        max_episodes = std::max(max_episodes, episode_statistics[iter].n_episodes);
    }

    for (int k=0; k<max_episodes; ++k) {
        int cnt = 0;
        real steps = 0;
        real total_reward = 0;
        real utility = 0;
        for (int iter=0; iter<n_iter; ++iter) {
            if (episode_statistics[iter].n_episodes > k){
                total_reward +=  episode_statistics[iter].total_reward[k];
                utility +=  episode_statistics[iter].utility[k];
                steps +=  episode_statistics[iter].steps[k];
                cnt++;
            }
        }
        real icnt = 1.0 / (real) cnt;
        printf ("%f %f %f %d #EPISODE\n",
                total_reward * icnt,
                utility * icnt,
                steps * icnt,
                cnt);
    }

    return 0;
}


/// This function views looks at a standard MDP-like thing.
bool EvaluateGeneral(Environment<Vector, int>& environment, 
                     real action_randomness,
                     FactoredPredictorRL<Vector, int>* factored_predictor,
                     real discount_factor,
                     EpisodeStatistics& episode_statistics,
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
    int steps = 0;
    real utility = 0;
    real total_reward = 0;
    real total_discount = 1;
    for (int t=0; t<T; ++t) {
        int action = last_action;
        if (!running) {
            factored_predictor->Observe(environment.getState());
            environment.Reset();
            factored_predictor->Reset();
            running = true;
            //printf("Terminated\n");
            if (t) {
                episode_statistics.AddEpisode(total_reward, utility, steps);
                steps = 0;
                utility = 0;
                total_reward = 0;
                total_discount = 1;
            }
        }

        real epsilon = action_randomness / (real) T;
        if (urandom() < epsilon) {
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
        //printf ("r:%f a:%d x:", reward, action); observation.print(stdout);
        real p = factored_predictor->Observe(action, observation, reward);
        //real td_error = factored_predictor->QLearning(0.1, discount_factor);
        real td_error = factored_predictor->Sarsa(0.1, 0.95, action_randomness);
        //assert(p==obs_probs[observation]);
        statistics.probability[t] += p;
        statistics.reward[t] += reward;
        statistics.error[t] += td_error;
        steps++;
        total_reward += reward;
        utility += total_discount * reward;
        total_discount *= discount_factor;
        //printf ("%d %d %f\n", action, observation, reward);
    }

    return true;
}

#endif
