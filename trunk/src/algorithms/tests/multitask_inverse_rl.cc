/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
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
#include "DiscreteChain.h"
#include "OptimisticTask.h"
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
#include "MWAL.h"
#include "PolicyBelief.h"
#include "RewardPolicyBelief.h"
#include "PolicyRewardBelief.h"
#include "PopulationPolicyRewardBelief.h"
#include "EasyClock.h"

#include <cstring>
#include <getopt.h>

struct EpisodeStatistics
{
    real total_reward;
    real discounted_reward;
    int steps;
    real mse;
    int n_runs;
    EpisodeStatistics()
        : total_reward(0.0),
          discounted_reward(0.0),
          steps(0),
          mse(0),
          n_runs(0)
    {

    }
};

#define N_COMPARISONS 4

struct Statistics
{
    std::vector<EpisodeStatistics> ep_stats;
    std::vector<real> reward;
    std::vector<int> n_runs;
    Vector DV;
    Matrix DV_total;
};


enum SamplerType {
    INVALID = 0x0,
    METROPOLIS,
    MONTE_CARLO
};

struct Sampler
{
    SamplerType type;
    int n_chains;
    int n_chain_samples;
};

/** Add a demonstration.

    \param n_steps maximum number of steps
    \param[out] demonstrations list of demonstrations
    \param epsilon-greedy/temperature parameter
    \param discount discount factor

    Returns the value iteration value function
*/
Vector AddDemonstration(uint n_steps,
                        Demonstrations<int, int>& demonstrations,
                        real epsilon, 
                        DiscreteEnvironment* environment, 
                        real discount)
{
    real accuracy = 1e-6;
    const DiscreteMDP* mdp = environment->getMDP(); 
    if (!mdp) {
        Serror("The environment must support the creation of an MDP\n");
        exit(-1);
    }

    ValueIteration value_iteration(mdp, discount);
    value_iteration.ComputeStateValues(accuracy);
    FixedSoftmaxPolicy softmax_policy(value_iteration.Q, 1.0 / epsilon);

    //std:: cout << "(value iteration)" << std::endl;
    
    environment->Reset();
    
    //std:: cout << "(running)" << std::endl;
    bool action_ok = true;

    for (uint step = 0; step < n_steps; ++step) {
        int state = environment->getState();
        int action = softmax_policy.SelectAction();

        action_ok = environment->Act(action);
        demonstrations.Observe(state, action);

        if (!action_ok) {
            break;
        }
    }
    demonstrations.NewEpisode();
    delete mdp;
    return value_iteration.V;
}

Statistics EvaluateAlgorithm (Demonstrations<int, int>& demonstrations,
                              std::vector<DiscreteEnvironment*>& environments,
                              real gamma,
                              real epsilon,
                              Sampler sampler,
                              real accuracy,
                              int iterations);

static const char* const environment_names_list = "{MountainCar, ContextBandit, RandomMDP, Gridworld, Chain, Optimistic}";

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
  --environment: {MountainCar, ContextBandit, RandomMDP, Gridworld, Chain, Optimistic}\n\
  --n_states:          number of states (usually there is no need to specify it)\n\
  --n_actions:         number of actions (usually there is no need to specify it)\n\
  --gamma:             reward discounting in [0,1]\n\
  --lambda:            eligibility trace parameter (for some algorithms)\n\
  --iterations:        number of iterations (for some algorithms)\n\
  --n_runs:            maximum number of runs\n\
  --n_demonstrations:        maximum number of episodes (ignored if < 0)\n\
  --randomness:        The higher this is, the more disperesed the distribution from which the tasks are drawn.\n\
  --n_tasks:           Specifies a finite number of tasks. For each run, n_tasks are sampled from the default task distribution. If not specified, or 0, the number of tasks is equal to the number of demonstrations.\n\
  --episode_steps:     maximum number of steps in each episode (ignored if <0)\n\
  --n_steps:           maximum number of total steps\n\
  --grid_size:         number of grid intervals for discretised environments\n\
  --maze_name:         (Gridworld) file name for the maze\n\
  --epsilon:           use epsilon-greedy with randomness in [0,1]\n\
  --n_chains:          number of chains [>0]\n\
  --n_chain_samples:  length of each chain [>0]\n\
  --Metropolis:        use Metropolis MCMC sampler\n\
  --MonteCarlo:        use plain MonteCarlo sampler\n\
  --accuracy:          estimation accuracy\n\
\n";


int main (int argc, char** argv) {

    int n_actions = 4;
    int n_states = 4;
    real gamma = 0.9;
    real lambda = 0.9;
    int iterations = 1000;
    real randomness = 0.1;
    real pit_value = -1.0;
    real goal_value = 1.0;
    real step_value = -0.01;
    real epsilon = 0.01;
    uint n_runs = 1000;
    uint n_steps = 100000;
    uint episode_steps = 1000;
    uint grid_size = 4;
    uint maze_width = 4;
    real transition_randomness = 0.01;
    uint n_tasks = 4;
    uint n_demonstrations = 100;
    Sampler sampler = {INVALID, 2, 10000};
    real accuracy = 1e-3;
    const char * environment_name = NULL;
    int max_samples = 4;
    char* maze_name = NULL;



    // options
    {
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
                {"n_demonstrations", required_argument, 0, 0}, //5
                {"n_steps", required_argument, 0, 0}, //6
                {"max_samples", required_argument, 0, 0}, //7
                {"multi-sample", no_argument, 0, 0}, //8
                {"maze_name", required_argument, 0, 0}, //9
                {"epsilon", required_argument, 0, 0}, //10
                {"maze_width", required_argument, 0, 0}, //11 - deprecated
                {"algorithm", required_argument, 0, 0}, //12
                {"environment", required_argument, 0, 0}, //13
                {"grid_size", required_argument, 0, 0}, //14
                {"randomness", required_argument, 0, 0}, //15
                {"episode_steps", required_argument, 0, 0}, //16
                {"iterations", required_argument, 0, 0}, //17
                {"n_chain_samples", required_argument, 0, 0}, //18
                {"n_chains", required_argument, 0, 0}, //19
                {"Metropolis", no_argument, 0, 0}, //20
                {"MonteCarlo", no_argument, 0, 0}, //21
                {"accuracy", required_argument, 0, 0}, //22
                {"n_tasks", required_argument, 0, 0}, //23
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
                case 5: n_demonstrations = atoi(optarg); break;
                case 6: n_steps = atoi(optarg); break;
                case 7: max_samples = atoi(optarg); break;
                case 8: printf("multi-sample not implented; ignored\n"); break;
                case 9: maze_name = optarg; break;
                case 10: epsilon = atof(optarg); break; 
                case 11: maze_width = atoi(optarg); break; // deprecated
                    //case 12: algorithm_name = optarg; break;
                case 13: environment_name = optarg; break;
                case 14: grid_size = atoi(optarg); break;
                case 15: randomness = atof(optarg); break;
                case 16: episode_steps = atoi(optarg); break;
                case 17: iterations = atoi(optarg); break;
                case 18: sampler.n_chain_samples = atoi(optarg); break;
                case 19: sampler.n_chains = atoi(optarg); break;
                case 20: sampler.type=METROPOLIS; break;
                case 21: sampler.type=MONTE_CARLO; break;
                case 22: accuracy = atof(optarg); assert(accuracy > 0); break;
                case 23: n_tasks = atoi(optarg); assert(n_tasks > 0); break;
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
    assert (n_demonstrations > 0);
    assert (n_steps > 0);
    assert (grid_size > 0);
    if (sampler.type == INVALID) {
        Serror ("invalid sampler specified\n");
        exit(-1);
    }
    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);



    
    RandomNumberGenerator* rng;
    
    //RandomNumberFile random_file("/home/olethros/projects/beliefbox/dat/r1e7.bin");
    //rng = (RandomNumberGenerator*) &random_file;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;
    rng->manualSeed(380731572);


    if (!environment_name) {
        std::cerr << "Please choose an environment from "
                  << environment_names_list << std::endl;
        exit(-1);
    }

    std::cout << "Starting test program" << std::endl;
    std::cout << "Starting evaluation" << std::endl;

    // }}} end load options


    // {{{ Set up statistics
    Statistics statistics;
    statistics.ep_stats.resize(n_demonstrations);
    statistics.reward.resize(n_steps);
    statistics.n_runs.resize(n_steps);
    for (uint i=0; i<n_steps; ++i) {
        statistics.reward[i] = 0;
        statistics.n_runs[i] = 0;
    }
    statistics.DV.Resize(N_COMPARISONS);
    statistics.DV.Clear();
    statistics.DV_total.Resize(n_runs, N_COMPARISONS);
    statistics.DV_total.Clear();
    // }}} end statistics set up



    for (uint run=0; run<n_runs; ++run) {
        logmsg("Run %d\n", 0);

        Demonstrations<int, int> demonstrations;
        std::vector<DiscreteEnvironment*> tasks(n_tasks);
        std::vector<DiscreteEnvironment*> environments(n_demonstrations);

        // sample the number of necessary environments
        for (uint i=0; i<n_tasks; ++i) {
            logmsg("Task %d\n", i);
            DiscreteEnvironment* environment = NULL;
            if (!strcmp(environment_name, "RandomMDP")) { 
                environment = new RandomMDP (n_actions,
                                             n_states,
                                             transition_randomness,
                                             step_value,
                                             pit_value,
                                             goal_value,
                                             rng,
                                             false);
            } else if (!strcmp(environment_name, "Gridworld")) {
                real pit_value = -1.0 * (1.0 - randomness) 
                    + urandom()* randomness;
                real step_value = -0.1 * (1.0 - randomness) 
                    + (urandom() - 0.5)* randomness;
                real goal_value = 1.0 * (1.0 - randomness)
                    + urandom()* randomness;
                environment = new Gridworld(maze_name,
                                            transition_randomness,
                                            pit_value,
                                            goal_value,
                                            step_value);
            } else if (!strcmp(environment_name, "ContextBandit")) { 
                environment = new ContextBandit(n_actions, 3, 4, rng);
            } else if (!strcmp(environment_name, "OneDMaze")) { 
                environment = new OneDMaze(n_states, rng);
            } else if (!strcmp(environment_name, "MountainCar")) { 
                MountainCar continuous_mountain_car;
                continuous_mountain_car.setRandomness(randomness);
                environment = new DiscretisedEnvironment<MountainCar> (continuous_mountain_car,  grid_size);
            } else if (!strcmp(environment_name, "Chain")) { 
                environment  = new DiscreteChain (n_states,
                                                  0.2,
                                                  0.2 * (1.0 - randomness) + urandom()*randomness,
                                                  1.0 * (1.0 - randomness) + urandom()*randomness);
            } else if (!strcmp(environment_name, "Optimistic")) { 
                environment = new OptimisticTask (0.1 * (1.0 - randomness) + urandom() * randomness, 0.01 * (1.0 - randomness) + urandom() * randomness);
            } else {
                fprintf(stderr, "Uknown environment %s\n", environment_name);
            }
            tasks[i] = environment;
            // making sure the number of states & actions is correct
            n_states = environment->getNStates();
            n_actions = environment->getNActions();
            
            std::cout <<  "# Creating environment: " << environment_name
                      << " with " << n_states << " states, "
                      << n_actions << " actions.\n";
        }

        // put at least one environment in each demonstration
        for (uint demonstration=0;
             demonstration<n_demonstrations;
             ++demonstration) {
            if (demonstration < n_tasks) {
                environments[demonstration] = tasks[demonstration];
            } else {
                int i = rng->discrete_uniform(n_tasks);
                environments[demonstration] = tasks[i];
            }
        }
        
        
        // {{{ Create data for demonstrations
        std::vector<Vector> demonstrator_value;
        logmsg("Creating data for demonstrations.\n");
        fflush(stdout);

        for (uint d = 0; d  < n_demonstrations; ++d) {
            DiscreteEnvironment* environment = environments[d];
            Vector V = AddDemonstration(episode_steps,
                                        demonstrations,
                                        epsilon,
                                        environment,
                                        gamma);
            demonstrator_value.push_back(V);
            //V.print(stdout);
        } // for d
        // }}} demo data

        logmsg("Evaluating.\n");
        fflush(stdout);

        // {{{ Evaluate methods
        Statistics run_statistics = EvaluateAlgorithm(demonstrations,
                                                      environments,
                                                      gamma,
                                                      epsilon,
                                                      sampler,
                                                      accuracy,
                                                      iterations);
        // }}}
        for (uint i=0; i<n_tasks; ++i) {
            delete tasks[i];
        }
    } // for runs
    
    return 0;
}
 
/*** Evaluate an algorithm

     episode_steps: maximum number of steps per episode. If negative, then ignore
     n_steps: maximun number of total steps. If negative, then ignore.
     n_demonstrations: number of demonstrations to do. Cannot be negative.
*/

Statistics EvaluateAlgorithm (Demonstrations<int, int>& demonstrations,
                              std::vector<DiscreteEnvironment*>& environments,
                              real gamma,
                              real epsilon,
                              Sampler sampler,
                              real accuracy,
                              int iterations)
{
    assert(accuracy >= 0);
    assert((int) demonstrations.size() == (int) environments.size());
    int n_demonstrations = environments.size();
    const DiscreteMDP* base_mdp = environments[0]->getMDP();
    
    if (!base_mdp) {
        Serror("The environment must support the creation of an MDP\n");
        exit(-1);
    }

    if (!base_mdp->Check()) {
        Serror("MDP model is nonsense\n");
        base_mdp->ShowModel();
    }

    int n_states = base_mdp->getNStates();
    int n_actions = base_mdp->getNActions();



    Statistics statistics;
    statistics.ep_stats.reserve(n_demonstrations); 

    // =================================== //
    // ||   T R A I N   M O D E L S     || //
    // =================================== //
    
    double start_time;

    // ------- Population model ------ //
    logmsg("Training PPRB\n"); fflush(stdout);
    start_time = GetCPU();
    PopulationPolicyRewardBelief pprb(1.0, gamma, *base_mdp);
    pprb.setAccuracy(accuracy);
    pprb.MonteCarloSampler(demonstrations, sampler.n_chain_samples);
    std::vector<FixedDiscretePolicy*> pprb_policies = pprb.getPolicies();
    printf("%f # T_PPRB\n", GetCPU() - start_time);

    // -------- MWAL -------- //
    logmsg("Training MWAL\n"); fflush(stdout);
    start_time = GetCPU();
    DiscreteMDP* computation_mdp = environments[0]->getMDP();
    MWAL mwal(n_states, n_actions, gamma);
    mwal.CalculateFeatureCounts(demonstrations);
    mwal.Compute(*computation_mdp, gamma, accuracy);//, iterations);
    printf("%f # T_MWAL\n", GetCPU() - start_time);
    delete computation_mdp;
    
    // ----- PRB: Policy | Reward Belief ------ //
    logmsg("Training SPRB\n"); fflush(stdout);
    start_time = GetCPU();
    PolicyRewardBelief sprb(1.0, gamma, *base_mdp);
    if (sampler.type == METROPOLIS) {
        logmsg("Metropolis sampler: %d %d\n", sampler.n_chain_samples, sampler.n_chains);
        sprb.MHSampler(demonstrations,
                      sampler.n_chain_samples,
                      sampler.n_chains);
    } else if (sampler.type == MONTE_CARLO) {
        logmsg("Monte Carlo sampler: %d %d\n", sampler.n_chain_samples, sampler.n_chains);
        sprb.MonteCarloSampler(demonstrations, sampler.n_chain_samples);
    } 
    printf("%f # T_SPRB\n", GetCPU() - start_time);
    DiscretePolicy* sprb_policy = sprb.getPolicy();
    
    // -------- imitator -------- //
    logmsg("Training IMIT\n"); fflush(stdout);
    start_time = GetCPU();
    printf("%f # T_IMIT\n", GetCPU() - start_time);
    FixedDiscretePolicy imitating_policy(n_states, n_actions,
                                         demonstrations);


    // =================================== //
    // ||    T E S T   M O D E L S      || //
    // =================================== //

    statistics.DV.Resize(N_COMPARISONS);
    for (int i=0; i<N_COMPARISONS; ++i) {
        statistics.DV(i) = 0;
    }
    for (int d = 0; d<n_demonstrations; ++ d) {
        DiscreteEnvironment* environment = environments[d];
        std:: cout << "# evaluating..." << environment->Name() << std::endl;
        const DiscreteMDP* mdp = environment->getMDP(); 
        ValueIteration value_iteration(mdp, gamma);
        value_iteration.ComputeStateValues(accuracy);        
        Vector& V_opt = value_iteration.V;
        for (int i=0; i<n_states; ++i) {
            printf ("%f ", V_opt(i));
        }
        printf (" # OPT\n");
        
        // -- IMIT -- //
        PolicyEvaluation imitating_evaluator(&imitating_policy, mdp, gamma);
        imitating_evaluator.ComputeStateValues(accuracy);
        for (int i=0; i<n_states; ++i) {
            printf ("%f ", imitating_evaluator.getValue(i));
        }   
        printf ("# IMIT \n");

        // -- MWAL -- //
        PolicyEvaluation mwal_evaluator(&mwal.mean_policy, mdp, gamma);
        mwal_evaluator.ComputeStateValues(accuracy);
        for (int i=0; i<n_states; ++i) {
            printf ("%f ", mwal_evaluator.getValue(i));
        }
        printf ("# MWAL\n");



        // -- SPRB -- //
        PolicyEvaluation sprb_evaluator(sprb_policy, mdp, gamma);
        sprb_evaluator.ComputeStateValues(accuracy);
        for (int i=0; i<n_states; ++i) {
            printf ("%f ", sprb_evaluator.getValue(i));
        }
        printf ("# SPRB \n");

        // -- PPRB -- //
        Vector V_PPRB(n_states);
        PolicyEvaluation pprb_evaluator(pprb_policies[d],mdp, gamma);
        pprb_evaluator.ComputeStateValues(accuracy);
        for (int i=0; i<n_states; ++i) {
            printf ("%f ", pprb_evaluator.getValue(i));
        }
        printf ("# PPRB \n");
        delete mdp;
        
        statistics.DV(0) += (V_opt - imitating_evaluator.V).L1Norm();
        statistics.DV(1) += (V_opt - mwal_evaluator.V).L1Norm();
        statistics.DV(2) += (V_opt - sprb_evaluator.V).L1Norm();
        statistics.DV(3) += (V_opt - pprb_evaluator.V).L1Norm();
        fflush(stdout);
        fflush(stderr);
    }
    statistics.DV /= (real) n_demonstrations;

#if 1
    printf ("%f %f %f %f # DV run\n", 
            statistics.DV(0),
            statistics.DV(1),
            statistics.DV(2),
            statistics.DV(3));
#endif        
    fflush(stdout);
    fflush(stderr);

    delete base_mdp;
    delete sprb_policy;
    for (int d=0; d<n_demonstrations; ++d) {
        delete pprb_policies[d];
    }
    return statistics;
}



#endif
