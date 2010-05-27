/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
// ------------ Environments -------------
// MDPs
#include "RandomMDP.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "ContextBandit.h"
// POMDPs
#include "OneDMaze.h"
#include "DiscretisedContinuousChain.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "DiscretisedEnvironment.h"

//-----------------------------------------

//------------- Algorithms ----------------
#include "PolicyEvaluation.h"
#include "PolicyIteration.h"
#include "ValueIteration.h"
#include "DiscretePolicy.h"
#include "Environment.h"
#include "ExplorationPolicy.h"
#include "Sarsa.h"
#include "QLearning.h"
#include "QLearningDirichlet.h"
#include "ModelBasedRL.h"
#include "ModelCollectionRL.h"
#include "ContextBanditGaussian.h"
#include "DiscreteMDPCollection.h"
#include "ContextBanditCollection.h"
//-------------------------------------------

#include "RandomNumberFile.h"
#include "MersenneTwister.h"
#include <cstring>

struct EpisodeStatistics
{
    real total_reward;
    real discounted_reward;
    int steps;
    real mse;
    EpisodeStatistics()
        : total_reward(0.0),
          discounted_reward(0.0),
          steps(0),
          mse(0)
    {

    }
};

struct Statistics
{
    std::vector<EpisodeStatistics> ep_stats;
    std::vector<real> reward;
};

Statistics EvaluateAlgorithm(uint n_steps,
                             uint n_episodes,
                             OnlineAlgorithm<int,int>* algorithm,
                             DiscreteEnvironment* environment,
                             real gamma);

int main (int argc, char** argv)
{
    int n_actions = 4;
    int n_original_states = 4;
    real gamma = 0.9;
    real lambda = 0.9;
    real alpha = 0.01;
    real randomness = 0.01;
    real pit_value = -1.0;
    real goal_value = 1.0;
    real step_value = -0.01;
    real epsilon = 0.01;
    uint n_runs = 1000;
    uint n_episodes = 1000;
    uint n_steps = 100;

    if (argc != 11) {
        for (int i=0; i<argc; ++i) {
            printf("arg %d: %s\n", i, argv[i]);
        }
        std::cerr << "Usage: online_algorithms n_states n_actions gamma lambda randomness n_runs n_episodes n_steps algorithm environment\n";
        return -1;
    }
    n_original_states = atoi(argv[1]);
    assert (n_original_states > 0);

    n_actions = atoi(argv[2]);
    assert (n_actions > 0);

    gamma = atof(argv[3]);
    assert (gamma >= 0 && gamma <= 1);

    lambda = atof(argv[4]);
    assert (lambda >= 0 && lambda <= 1);

    randomness = atof(argv[5]);
    assert (randomness >= 0 && randomness <= 1);
    
    n_runs = atoi(argv[6]);
    assert (n_runs > 0);

    n_episodes = atoi(argv[7]);
    assert (n_episodes > 0);

    n_steps = atoi(argv[8]);
    assert (n_steps > 0);

    char* algorithm_name = argv[9];
    char* environment_name = argv[10];

    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);

    DiscreteMDPCounts* discrete_mdp = NULL;
    
    RandomNumberGenerator* rng;
    
    //RandomNumberFile random_file("/home/olethros/projects/beliefbox/dat/r1e7.bin");
    //rng = (RandomNumberGenerator*) &random_file;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;

    std::cout << "Starting test program" << std::endl;
    
    std::cout << "Starting evaluation" << std::endl;
    // remember to use n_runs
    Statistics statistics;
    statistics.ep_stats.resize(n_episodes);
    //statistics.reward.resize(n_steps);
    for (uint i=0; i<statistics.ep_stats.size(); ++i) {
        statistics.ep_stats[i].total_reward = 0.0;
        statistics.ep_stats[i].discounted_reward = 0.0;
        statistics.ep_stats[i].steps = 0;
        statistics.ep_stats[i].mse = 0;
    }
    for (uint run=0; run<n_runs; ++run) {
        std::cout << "Run: " << run << " - Creating environment.." << std::endl;
        std::cerr << "Run: " << run << " - Creating environment.." << std::endl;
        DiscreteEnvironment* environment = NULL;

        RandomMDP* random_mdp = NULL;
        OneDMaze* one_d_maze = NULL; 
        DiscretisedContinuousChain* discretised_chain = NULL; 
        MountainCar continuous_mountain_car;
        DiscretisedEnvironment<MountainCar>* mountain_car = NULL;
        Pendulum continuous_pendulum;
        DiscretisedEnvironment<Pendulum>* pendulum = NULL;
        Gridworld* gridworld = NULL;
        ContextBandit* context_bandit = NULL;
        if (!strcmp(environment_name, "RandomMDP")) { 
            random_mdp = new RandomMDP (n_actions,
                                        n_original_states,
                                        randomness,
                                        step_value,
                                        pit_value,
                                        goal_value,
                                        rng,
                                        false);
            environment = random_mdp;
        } else if (!strcmp(environment_name, "Gridworld")) { 
            gridworld = new Gridworld("/home/olethros/projects/beliefbox/dat/maze2",  8, 8);
            environment = gridworld;
        } else if (!strcmp(environment_name, "ContextBandit")) { 
            context_bandit = new ContextBandit(n_actions, 3, 4, rng);
            environment = context_bandit;
        } else if (!strcmp(environment_name, "OneDMaze")) { 
            one_d_maze = new OneDMaze(n_original_states, rng);
            environment = one_d_maze;
        } else if (!strcmp(environment_name, "DiscretisedContinuousChain")) { 
            discretised_chain = new DiscretisedContinuousChain(n_original_states);
            environment = discretised_chain;
        } else if (!strcmp(environment_name, "MountainCar")) { 
            mountain_car = new DiscretisedEnvironment<MountainCar>(continuous_mountain_car, n_original_states);
            environment = mountain_car;
        } else if (!strcmp(environment_name, "Pendulum")) { 
            pendulum = new DiscretisedEnvironment<Pendulum>(continuous_pendulum, n_original_states);
            environment = pendulum;
        } else {
            fprintf(stderr, "Uknown environment %s\n", environment_name);
        }


        // making sure the number of states & actions is correct
        int n_states = environment->getNStates();
        n_actions = environment->getNActions();
        
        std::cout <<  "Creating environment: " << environment_name
                  << " with " << n_states << "states, "
                  << n_actions << " actions.\n";

        //std::cout << "Creating exploration policy" << std::endl;
        VFExplorationPolicy* exploration_policy = NULL;
        exploration_policy = new EpsilonGreedy(n_actions, epsilon);
    
    
        //std::cout << "Creating online algorithm" << std::endl;
        OnlineAlgorithm<int, int>* algorithm = NULL;
        MDPModel* model = NULL;
        //Gridworld* g2 = gridworld;
        if (!strcmp(algorithm_name, "Sarsa")) { 
            algorithm = new Sarsa(n_states,
                                  n_actions,
                                  gamma,
                                  lambda,
                                  alpha,
                                  exploration_policy);
        } else if (!strcmp(algorithm_name, "QLearning")) { 
            algorithm = new QLearning(n_states,
                                      n_actions,
                                      gamma,
                                      lambda,
                                      alpha,
                                      exploration_policy);
        } else if (!strcmp(algorithm_name, "QLearningDirichlet")) { 
            algorithm = new QLearningDirichlet(n_states,
                                               n_actions,
                                               gamma,
                                               lambda,
                                               alpha,
                                               exploration_policy);
        } else if (!strcmp(algorithm_name, "Model")) {
            discrete_mdp =  new DiscreteMDPCounts(n_states, n_actions);
            model= (MDPModel*) discrete_mdp;

            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         model);
        } else if (!strcmp(algorithm_name, "ContextBanditGaussian")) {
            model= (MDPModel*)
                new ContextBanditGaussian(n_states,
                                          n_actions,
                                          0.5, 0.0, 1.0);
            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         model,
                                         false);
        } else if (!strcmp(algorithm_name, "Aggregate")) {
            model= (MDPModel*)
                new ContextBanditAggregate(false, 3, 2,
                                           n_states, 4,
                                           n_actions,
                                           0.5, 0.0, 1.0);
            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         model,
                                         false);
        } else if (!strcmp(algorithm_name, "Collection")) {
            DiscreteMDPCollection* collection = 
                new DiscreteMDPCollection(*gridworld, 
                                          8,
                                          n_states,
                                          n_actions);
            model= (MDPModel*) collection;

            algorithm = new ModelCollectionRL(n_states,
                                              n_actions,
                                              gamma,
                                              epsilon,
                                              collection,
                                              true);
        } else if (!strcmp(algorithm_name, "ContextBanditCollection")) {
            ContextBanditCollection* collection = 
                new ContextBanditCollection(8,
                                            n_states,
                                            n_actions,
                                            0.5, 0.0, 1.0);
            model= (MDPModel*) collection;

            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         collection,
                                         false);
        } else {
            Serror("Unknown algorithm: %s\n", algorithm_name);
        }

        
        //std::cerr << "run : " << run << std::endl;
        Statistics run_statistics = EvaluateAlgorithm(n_steps, n_episodes, algorithm, environment, gamma);
        for (uint i=0; i<statistics.ep_stats.size(); ++i) {
            statistics.ep_stats[i].total_reward += run_statistics.ep_stats[i].total_reward;
            statistics.ep_stats[i].discounted_reward += run_statistics.ep_stats[i].discounted_reward;
            statistics.ep_stats[i].steps += run_statistics.ep_stats[i].steps;
            statistics.ep_stats[i].mse += run_statistics.ep_stats[i].mse;
        }

        //for (uint i=0; i<statistics.reward.size(); ++i) {
        //statistics.reward[i] += run_statistics.reward[i];
        //}
        if (model) {
#if 0
            if (discrete_mdp) {
                real threshold = 10e-6; //0;
                int max_iter = 10;//100;
                DiscreteMDP* mean_mdp = discrete_mdp->getMeanMDP();
                PolicyIteration MPI(mean_mdp, gamma);
                ValueIteration MVI(mean_mdp, gamma);
                MPI.ComputeStateValues(threshold, max_iter);
                MVI.ComputeStateValues(threshold, max_iter);
                FixedDiscretePolicy* policy = MVI.getPolicy();
                PolicyEvaluation MPE(policy, mean_mdp, gamma);
                MPE.ComputeStateValues(threshold, max_iter);

                Vector hV(n_states);
                Vector hU(n_states);
                int n_samples = 1;
                Vector Delta(n_samples);
                for (int i=0; i<n_samples; ++i) {
                    DiscreteMDP* sample_mdp = discrete_mdp->generate();
                    PolicyEvaluation PE(policy, sample_mdp, gamma);
                    PE.ComputeStateValues(threshold, max_iter);
                    ValueIteration VI(sample_mdp, gamma);
                    VI.ComputeStateValues(threshold, max_iter);
                    //VI.ComputeStateActionValues(threshold, max_iter);
                    Delta[i] = 0.0;
                    for (int s=0; s<n_states; ++s) {
                        hV[s] += PE.getValue(s);
                        hU[s] += VI.getValue(s);
                        Delta[i] += fabs(MVI.getValue(s) - hV[s] / (real) i);
                    }
                    printf ("%f #delta\n", Delta[i]);
                    delete sample_mdp;
                    
                }
                real inv_n = 1.0 / (real) n_samples;
                for (int s=0; s<n_states; ++s) {
                    hV[s] *= inv_n;
                    hU[s] *= inv_n;
                    printf ("V[%d] = %f %f %f | %f %f\n",
                    s, MPI.getValue(s), MPE.getValue(s), MVI.getValue(s), hV[s], hU[s]);
                    //printf ("%f %f #hV\n", MVI.getValue(s) hU[s]);

                }
                delete mean_mdp;
                delete policy;

            }
#endif
            delete model;
            model = NULL;
        }
        //delete environment;
        delete gridworld;
        delete one_d_maze;
        delete discretised_chain;
        delete random_mdp;
        delete context_bandit;
        delete algorithm;
        delete exploration_policy;
    }
    

    for (uint i=0; i<statistics.ep_stats.size(); ++i) {
        statistics.ep_stats[i].total_reward /= (float) n_runs;
        statistics.ep_stats[i].discounted_reward /= (float) n_runs;
        statistics.ep_stats[i].steps /= n_runs;
        statistics.ep_stats[i].mse /= n_runs;
        std::cout << i << "i" << std::endl;

        std::cout << statistics.ep_stats[i].total_reward << " "
                  << statistics.ep_stats[i].discounted_reward << "# REWARD"
                  << std::endl;
        std::cout << statistics.ep_stats[i].steps << " "
                  << statistics.ep_stats[i].mse << "# MSE"
                  << std::endl;
    }

    //for (uint i=0; i<statistics.reward.size(); ++i) {
    //        statistics.reward[i] /= (float) n_runs;
    //        std::cout << statistics.reward[i] << " # INST_PAYOFF"
    //                  << std::endl;
    //    }

    std::cout << "Done" << std::endl;


    
    return 0;
}

Statistics EvaluateAlgorithm (uint n_steps,
                             uint n_episodes,
                             OnlineAlgorithm<int, int>* algorithm,
                             DiscreteEnvironment* environment,
                             real gamma)
{
    std:: cout << "evaluating..." << environment->Name() << std::endl;
    
    std:: cout << "(value iteration)" << std::endl;

    //int n_states = environment->getNStates();
    //int n_actions = environment->getNActions();

    Statistics statistics;
    statistics.ep_stats.resize(n_episodes);
    //statistics.reward.resize(n_episodes*n_steps);

    real discount = 1.0;
    int current_time = 0;
    environment->Reset();
    std::cout << "(running)" << std::endl;
    int spin_state = 0;
    //const char* spinner ="|/-\\";
    const char* spinner ="_,.-'`'-.,_";
    int max_spin_state = strlen(spinner);
    for (uint episode = 0; episode < n_episodes; ++episode) {
        fprintf(stderr, "\r%c%d", spinner[spin_state], 100*episode/n_episodes);
        spin_state = (spin_state + 1) % max_spin_state;
        statistics.ep_stats[episode].total_reward = 0.0;
        statistics.ep_stats[episode].discounted_reward = 0.0;
        statistics.ep_stats[episode].steps = 0;
        discount = 1.0;
        environment->Reset();
        algorithm->Reset();
        uint t;
        for (t=0; t < n_steps; ++t) {
            int state = environment->getState();
            real reward = environment->getReward();
            //statistics.reward[current_time] = reward;
            statistics.ep_stats[episode].total_reward += reward;
            statistics.ep_stats[episode].discounted_reward += discount * reward;
            discount *= gamma;
            //std::cout << "Acting!\n";
            int action = algorithm->Act(reward, state);
            //std::cout << "s:" << state << " r:" << reward << " a:" << action << std::endl;
            bool action_ok = environment->Act(action);
            if (!action_ok) {
                break;
            }
            current_time++;
        }
        statistics.ep_stats[episode].steps += t;
        statistics.ep_stats[episode].mse  = 0;
    }
    fprintf(stderr, "\n");
    std::cout << "(done)" << std::endl;
    //std::cout << "REAL MODEL\n";
    //mdp->ShowModel();
    
    return statistics;
}

#endif
