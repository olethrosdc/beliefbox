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
#include "PolicyEvaluation.h"
#include "ValueIteration.h"
#include "RandomMDP.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"
#include "Environment.h"
#include "ExplorationPolicy.h"
#include "Sarsa.h"
#include "QLearning.h"
#include "QLearningDirichlet.h"
#include "ModelBasedRL.h"
#include "ContextBanditGaussian.h"
#include "ContextBandit.h"
#include "DiscreteMDPCollection.h"
#include "ContextBanditCollection.h"
#include "RandomNumberFile.h"

struct EpisodeStatistics
{
    real total_reward;
    real discounted_reward;
    int steps;
    real mse;
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
    int n_states = 4;
    real gamma = 0.9;
    real lambda = 0.9;
    real alpha = 0.1;
    real randomness = 0.1;
    real pit_value = -1.0;
    real goal_value = 1.0;
    real step_value = -0.01;
    real epsilon = 0.1;
    uint n_runs = 1000;
    uint n_episodes = 1000;
    uint n_steps = 100;

    if (argc != 10) {
        std::cerr << "Usage: online_algorithms n_states n_actions gamma lambda randomness n_runs n_episodes n_steps algorithm\n";
        return -1;
    }
    n_states = atoi(argv[1]);
    assert (n_states > 0);

    n_actions = atoi(argv[2]);
    assert (n_actions > 0);

    gamma = atof(argv[3]);
    assert (gamma > 0 && gamma <= 1);

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

    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);
   
    RandomNumberFile random_file("/home/olethros/projects/beliefbox/dat/r1e7.bin");
    RandomNumberGenerator* rng = (RandomNumberGenerator*) &random_file;

    std::cout << "Starting test program" << std::endl;
    
    std::cout << "Starting evaluation" << std::endl;
    // remember to use n_runs
    Statistics statistics;
    statistics.ep_stats.resize(n_episodes);
    statistics.reward.resize(n_episodes*n_steps);
    for (uint i=0; i<statistics.ep_stats.size(); ++i) {
        statistics.ep_stats[i].total_reward = 0.0;
        statistics.ep_stats[i].discounted_reward = 0.0;
        statistics.ep_stats[i].steps = 0;
        statistics.ep_stats[i].mse = 0;
    }
    for (uint run=0; run<n_runs; ++run) {
        std::cout << "Run: " << run << " - Creating environment.." << std::endl;
        DiscreteEnvironment* environment = NULL;
#if 0
        environment = new RandomMDP (n_actions,
                                     n_states,
                                     randomness,
                                     step_value,
                                     pit_value,
                                     goal_value,
                                     rng,
                                     false);
#elseif 0
        Gridworld* gridworld= new Gridworld("/home/olethros/projects/beliefbox/dat/maze3",
                                            16, 16);
        environment = gridworld;
#else
        ContextBandit* context_bandit = new ContextBandit(n_actions, 3, 8, rng);
        environment = context_bandit;
#endif
    
        // making sure the number of states & actions is correct
        n_states = environment->getMDP()->GetNStates();
        n_actions = environment->getMDP()->GetNActions();


        //std::cout << "Creating exploration policy" << std::endl;
        VFExplorationPolicy* exploration_policy = NULL;
        exploration_policy = new EpsilonGreedy(n_actions, epsilon);
    
    
        //std::cout << "Creating online algorithm" << std::endl;
        OnlineAlgorithm<int, int>* algorithm = NULL;
        MDPModel* model = NULL;
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
            model= (MDPModel*)
                new DiscreteMDPCounts(n_states,
                                      n_actions);
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
        } else if (!strcmp(algorithm_name, "Collection")) {
#if 1
            //new DiscreteMDPCollection(8,
            model= (MDPModel*) 
                new ContextBanditCollection(8,
                                          n_states,
                                          n_actions,
                                          0.5, 0.0, 1.0);
#else
            model= (MDPModel*)
                new DiscreteMDPCollection(*gridworld, 
                                          4,
                                          n_states,
                                          n_actions);
#endif
            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         model,
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
        for (uint i=0; i<statistics.reward.size(); ++i) {
            statistics.reward[i] += run_statistics.reward[i];
        }
        if (model) {
            delete model;
            model = NULL;
        }
        delete environment;
        delete algorithm;
        delete exploration_policy;
    }
    

    for (uint i=0; i<statistics.ep_stats.size(); ++i) {
        statistics.ep_stats[i].total_reward /= (float) n_runs;
        statistics.ep_stats[i].discounted_reward /= (float) n_runs;
        statistics.ep_stats[i].steps /= n_runs;
        statistics.ep_stats[i].mse /= n_runs;
        std::cout << statistics.ep_stats[i].total_reward << " "
                  << statistics.ep_stats[i].discounted_reward << "# REWARD"
                  << std::endl;
        std::cout << statistics.ep_stats[i].steps << " "
                  << statistics.ep_stats[i].mse << "# MSE"
                  << std::endl;
    }

    for (uint i=0; i<statistics.reward.size(); ++i) {
        statistics.reward[i] /= (float) n_runs;
        std::cout << statistics.reward[i] << " # INST_PAYOFF"
                  << std::endl;
    }

    std::cout << "Done" << std::endl;


    
    return 0;
}

Statistics EvaluateAlgorithm (uint n_steps,
                             uint n_episodes,
                             OnlineAlgorithm<int, int>* algorithm,
                             DiscreteEnvironment* environment,
                             real gamma)
{
    std:: cout << "evaluating..." << std::endl;
 
    const DiscreteMDP* mdp = environment->getMDP(); 
    ValueIteration value_iteration(mdp, gamma);
    if (!mdp) {
        Serror("The environment must support the creation of an MDP\n");
        exit(-1);
    }
    std:: cout << "(value iteration)" << std::endl;
    value_iteration.ComputeStateActionValues(10e-6,10);
    int n_states = mdp->GetNStates();
    int n_actions = mdp->GetNActions();

    Statistics statistics;
    statistics.ep_stats.resize(n_episodes);
    statistics.reward.resize(n_episodes*n_steps);

    real discount = 1.0;
    int current_time = 0;
    environment->Reset();
    std:: cout << "(running)" << std::endl;
    for (uint episode = 0; episode < n_episodes; ++episode) {
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
            statistics.reward[current_time] = reward;
            statistics.ep_stats[episode].total_reward += reward;
            statistics.ep_stats[episode].discounted_reward += discount * reward;
            discount *= gamma;
            int action = algorithm->Act(reward, state);
            //std::cout << "s:" << state << " r:" << reward << " a:" << action << std::endl;
            bool action_ok = environment->Act(action);
            if (!action_ok) {
                break;
            }
            current_time++;
        }
        statistics.ep_stats[episode].steps += t;
        real sse = 0.0;
        for (int i=0; i<n_states; i++) {
            for (int a=0; a<n_actions; a++) {
                real V =  value_iteration.getValue(i, a);
                real hV = algorithm->getValue(i, a);
                //printf ("Q(%d, %d) = %f ~ %f\n", i, a, V, hV);
                real err = V - hV;
                sse += err*err;
            }
        }
        statistics.ep_stats[episode].mse += sse /((real) (n_states*n_actions));
    }

    //std::cout << "REAL MODEL\n";
    //mdp->ShowModel();
    
    return statistics;
}

#endif
