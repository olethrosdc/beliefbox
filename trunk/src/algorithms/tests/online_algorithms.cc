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
#include "PolicyIteration.h"
#include "ValueIteration.h"
#include "RandomMDP.h"
#include "Gridworld.h"
#include "OneDMaze.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"
#include "Environment.h"
#include "ExplorationPolicy.h"
#include "Sarsa.h"
#include "QLearning.h"
#include "QLearningDirichlet.h"
#include "ModelBasedRL.h"
#include "SampleBasedRL.h"
#include "ModelCollectionRL.h"
#include "ContextBanditGaussian.h"
#include "ContextBandit.h"
#include "DiscreteMDPCollection.h"
#include "ContextBanditCollection.h"
#include "RandomNumberFile.h"
#include "MersenneTwister.h"
#include "MountainCar.h"
#include "DiscretisedEnvironment.h"
#include "HQLearning.h"
#include <cstring>
#include <getopt.h>

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
static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
    --algorithm:   {QLearning, Model, Sarsa, Sampling}\n\
    --environment: {MountainCar, ContextBandit, RandomMDP, Gridworld}\n\
    --n_states:    number of states (usually there is no need to specify it)\n\
    --n_actions:   number of actions (usually there is no need to specify it)\n\
    --gamma:       reward discounting in [0,1]\n\
    --lambda:      eligibility trace parameter (for some algorithms)\n\
    --randomness:  environment randomness\n\
    --n_runs:      maximum number of runs\n\
    --n_episodes:  maximum number of episodes\n\
    --n_steps:     maximum number of steps in each episode\n\
    --grid_size:   number of grid intervals for discretised environments\n\
    --maze_name:   (Gridworld) file name for the maze\n\
\n";


int main (int argc, char** argv)
{
    int n_actions = 4;
    int n_states = 4;
    real gamma = 0.9;
    real lambda = 0.9;
    real alpha = 0.1;
    real randomness = 0.01;
    real pit_value = -1.0;
    real goal_value = 1.0;
    real step_value = -0.01;
    real epsilon = 0.01;
    uint n_runs = 1000;
    uint n_episodes = 1000;
    uint n_steps = 100;
    uint grid_size = 4;
    uint maze_height = 4;
    uint maze_width = 4;
    const char * algorithm_name = "QLearning";
    const char * environment_name = "Gridworld";


    int max_samples = 4;
    char* maze_name = NULL;
    {
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"n_states", required_argument, 0, 0}, //0
                {"n_actions", required_argument, 0, 0}, //1
                {"gamma", required_argument, 0, 0}, //2
                {"lambda", required_argument, 0, 0}, //3
                {"n_runs", required_argument, 0, 0}, //4
                {"n_episodes", required_argument, 0, 0}, //5
                {"n_steps", required_argument, 0, 0}, //6
                {"max_samples", required_argument, 0, 0}, //7
                {"multi-sample", no_argument, 0, 0}, //8
                {"maze_name", required_argument, 0, 0}, //9
                {"maze_height", required_argument, 0, 0}, //10 - deprecataed
                {"maze_width", required_argument, 0, 0}, //11 - deprecated
                {"algorithm", required_argument, 0, 0}, //12
                {"environment", required_argument, 0, 0}, //13
                {"grid_size", required_argument, 0, 0}, //14
                {"randomness", required_argument, 0, 0}, //15
                {0, 0, 0, 0}
            };
            c = getopt_long (argc, argv, "",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
            case 0:
#if 0
                printf ("option %s (%d)", long_options[option_index].name, option_index);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
#endif
                switch (option_index) {
                case 0: n_states = atoi(optarg); break;
                case 1: n_actions = atoi(optarg); break;
                case 2: gamma = atof(optarg); break;
                case 3: lambda = atof(optarg); break;
                case 4: n_runs = atoi(optarg); break;
                case 5: n_episodes = atoi(optarg); break;
                case 6: n_steps = atoi(optarg); break;
                case 7: max_samples = atoi(optarg); break;
                case 8: printf("multi-sample not implented; ignored\n"); break;
                case 9: maze_name = optarg; break;
                case 10: maze_height = atoi(optarg); break; // deprecated
                case 11: maze_width = atoi(optarg); break; // deprecated
                case 12: algorithm_name = optarg; break;
                case 13: environment_name = optarg; break;
                case 14: grid_size = atoi(optarg); break;
                case 15: randomness = atof(optarg); break;
                default:
                    fprintf (stderr, "%s", help_text);
                    exit(0);
                    break;
                }
                break;
            case '0':
            case '1':
            case '2':
                if (digit_optind != 0 && digit_optind != this_option_optind)
                    printf ("digits occur in two different argv-elements.\n");
            digit_optind = this_option_optind;
            printf ("option %c\n", c);
            break;
            default:
                std::cout << help_text;
                exit (-1);
            }
        }
	
        if (optind < argc) {
            printf ("non-option ARGV-elements: ");
            while (optind < argc) {
                printf ("%s ", argv[optind++]);
                
            }
            printf ("\n");
        }
    }

    assert (n_states > 0);
    assert (n_actions > 0);
    assert (gamma >= 0 && gamma <= 1);
    assert (lambda >= 0 && lambda <= 1);
    assert (randomness >= 0 && randomness <= 1);
    assert (n_runs > 0);
    assert (n_episodes > 0);
    assert (n_steps > 0);
    assert (grid_size > 0);

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
        RandomMDP* random_mdp = new RandomMDP (n_actions,
                                               n_states,
                                               randomness,
                                               step_value,
                                               pit_value,
                                               goal_value,
                                               rng,
                                               false);
        OneDMaze* one_d_maze = new OneDMaze(n_states, rng);
        Gridworld* gridworld= new Gridworld(maze_name, randomness);
        ContextBandit* context_bandit = new ContextBandit(n_actions, 3, 4, rng);
        MountainCar continuous_mountain_car;
        DiscretisedEnvironment<MountainCar>* mountain_car
            = new DiscretisedEnvironment<MountainCar> (continuous_mountain_car,
                                                       n_states);
        if (!strcmp(environment_name, "RandomMDP")) { 
            environment = random_mdp;
        } else if (!strcmp(environment_name, "Gridworld")) { 
            environment = gridworld;
        } else if (!strcmp(environment_name, "ContextBandit")) { 
            environment = context_bandit;
        } else if (!strcmp(environment_name, "OneDMaze")) { 
            environment = one_d_maze;
        } else if (!strcmp(environment_name, "MountainCar")) { 
            environment = mountain_car;
        } else {
            fprintf(stderr, "Uknown environment %s\n", environment_name);
        }


        // making sure the number of states & actions is correct
        n_states = environment->getNStates();
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
        } else if (!strcmp(algorithm_name, "HQLearning")) { 
            algorithm = new HQLearning(
									   4,
									   n_states,
									   n_actions,
									   gamma,
									   lambda,
									   alpha,
									   0.01,
									   1.0);
        } else if (!strcmp(algorithm_name, "QLearningDirichlet")) { 
            algorithm = new QLearningDirichlet(n_states,
                                               n_actions,
                                               gamma,
                                               lambda,
                                               alpha,
                                               exploration_policy);
        } else if (!strcmp(algorithm_name, "Model")) {
            discrete_mdp =  new DiscreteMDPCounts(n_states, n_actions, 1.0 / (real) n_states);
            model= (MDPModel*) discrete_mdp;
            algorithm = new ModelBasedRL(n_states,
                                         n_actions,
                                         gamma,
                                         epsilon,
                                         model);
        } else if (!strcmp(algorithm_name, "Sampling")) {
            discrete_mdp =  new DiscreteMDPCounts(n_states, n_actions, 1.0 / (real) n_states);
            model= (MDPModel*) discrete_mdp;
            algorithm = new SampleBasedRL(n_states,
                                          n_actions,
                                          gamma,
                                          epsilon,
                                          model,
                                          max_samples);
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
            DiscreteMDPCollection* collection = NULL;
            if (environment == gridworld) {
                collection = new DiscreteMDPCollection(*gridworld, 
                                                       4,
                                                       n_states,
                                                       n_actions);
            } else {
                collection =  new DiscreteMDPCollection(2,
                                                        n_states,
                                                        n_actions);
            }
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
        for (uint i=0; i<statistics.reward.size(); ++i) {
            statistics.reward[i] += run_statistics.reward[i];
        }
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
        delete mountain_car;
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
    std:: cout << "evaluating..." << environment->Name() << std::endl;
    
#if 0
    const DiscreteMDP* mdp = environment->getMDP(); 
    //ValueIteration value_iteration(mdp, gamma);
    if (!mdp) {
        Serror("The environment must support the creation of an MDP\n");
        exit(-1);
    }
    std:: cout << "(value iteration)" << std::endl;
    //value_iteration.ComputeStateActionValues(10e-6,1000);
    int n_states = mdp->GetNStates();
    int n_actions = mdp->GetNActions();
#endif

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
        bool action_ok = true;
        uint t;
        for (t=0; t < n_steps; ++t) {
            int state = environment->getState();
            real reward = environment->getReward();
            statistics.reward[current_time] = reward;
            statistics.ep_stats[episode].total_reward += reward;
            statistics.ep_stats[episode].discounted_reward += discount * reward;
            discount *= gamma;
            //std::cout << "Acting!\n";
            int action = algorithm->Act(reward, state);
            //std::cout << "s:" << state << " r:" << reward << " a:" << action << std::endl;
            if (!action_ok) {
                break;
            }
            action_ok = environment->Act(action);

            current_time++;
        }
        statistics.ep_stats[episode].steps += t;
#if 0
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
#endif
    }

    //std::cout << "REAL MODEL\n";
    //mdp->ShowModel();
    
    return statistics;
}

#endif
