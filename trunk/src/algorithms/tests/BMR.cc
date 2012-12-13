/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "Grid.h"
#include "BasisSet.h"
#include "Random.h"
#include "MersenneTwister.h"
#include "Vector.h"
#include "Matrix.h"

// -- Discrete environments -- //
#include "RandomMDP.h"
#include "Gridworld.h"
#include "OneDMaze.h"
#include "DiscreteChain.h"
#include "RiverSwim.h"
#include "OptimisticTask.h"
#include "InventoryManagement.h"
#include "DoubleLoop.h"

// -- Continuous environments -- //
#include "MountainCar.h"
#include "Pendulum.h"
#include "PuddleWorld.h"
#include "DiscretisedEnvironment.h"

// -- Randomness -- //
#include "RandomNumberFile.h"
#include "MersenneTwister.h"

// -- Algorithms and models --//
#include "Grid.h"
#include "BasisSet.h"
#include "LinearModel.h"
#include "OnlineAlgorithm.h"
#include "BayesianMultivariateRegression.h"
#include "RepresentativeStateModel.h"
#include "BayesianMultivariate.h"


// -- Usual includes -- //
#include <cstring>
#include <getopt.h>

struct EpisodeStatistics
{
	real total_reward;
	real discounted_reward;
	int steps;
	real mse;
	int n_runs;
	EpisodeStatistics() : total_reward(0.0), 
						  discounted_reward(0.0),
						  steps(0),
						  mse(0),
						  n_runs(0)
	{}
};

struct Statistics
{
	std::vector<EpisodeStatistics> ep_stats;
	std::vector<real> reward;
	std::vector<int> n_runs;
};

void TrainingAlgorithm(real gamma,
					   int n_episodes,
					   BayesianMultivariate* algorithm, 
					   ContinuousStateEnvironment* environment,
					   RepresentativeStateModel<LinearModel<Vector,int>, Vector, int> *R);

Statistics EvaluateAlgorithm(int episode_steps,
							 int n_episodes,
							 uint n_steps,
							 BayesianMultivariate* algorithm,
							 ContinuousStateEnvironment* environment,
							 real gamma);

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
--environment:		{PuddleWorld, MountainCar, Pendulum}\n\
--n_states:			number of states (usually there is no need to specify it)\n\
--n_actions:		number of actions (usually there is no need to specify it)\n\
--gamma:			reward discounting in [0,1] (* 0.99)\n\
--randomness:		environment randomness (* 0.0)\n\
--n_runs:			maximum number of runs (* 1)\n\
--n_episodes:		maximum number of episodes (ignored if < 0)\n\
--episode_steps:    maximum number of steps in each episode (ignored if <0)\n\
--n_steps:			maximum number of total steps\n\
--b_functions:		use of basis functions or not (* 0)\n\
--grids:			number of grid intervals for discretised environments (*3)\n\
--epsilon:			use epsilon-greedy with randomness in [0,1] (* 0.01)\n\
--a:				first linear model parameter (*0.1) \n\
--N0:				second linear model parameter (*0.1) \n\
--n_samples:		number of collected samples (* 1000) \n\
\n\
* denotes default parameter\n\
\n";


int main(int argc, char* argv[])
{
	int n_actions			= 2;
	int n_states			= 5;
	real gamma				= 0.99;
	real randomness			= 0.0;
	real epsilon			= 0.1;
	real N0					= 0.1;
	real a					= 0.1;
	uint n_runs				= 1;
    int n_episodes			= 10;
	uint episode_steps		= 1000;
    uint n_steps			= 100000;
	uint b_functions		= 0;
    uint grids				= 3;
	int n_samples			= 100;
	int sampling			= 1;
	//Value iteration threshold
	real threshold = 0.00001;


	const char * environment_name = "MountainCar";
	{
		//options
		int c;
		int digit_optind = 0;
		while(1) {
			int this_option_optind = optind ? optind : 1;
			int option_index = 0;
			static struct option long_options[] = {
				{"n_states", required_argument, 0, 0},			//0
				{"n_actions", required_argument, 0, 0},			//1
				{"gamma", required_argument, 0, 0},				//2
				{"n_runs", required_argument, 0, 0},			//4
                {"n_episodes", required_argument, 0, 0},		//5
                {"n_steps", required_argument, 0, 0},			//6
                {"epsilon", required_argument, 0, 0},			//7
				{"N0", required_argument, 0, 0},				//8
				{"a", required_argument, 0, 0},					//9
                {"environment", required_argument, 0, 0},		//10
				{"b_functions", required_argument, 0, 0},		//11
                {"grids", required_argument, 0, 0},				//12
                {"randomness", required_argument, 0, 0},		//13
                {"episode_steps", required_argument, 0, 0},		//14
				{"n_samples", required_argument, 0, 0},			//15
				{"sampling", required_argument, 0, 0},			//16
				{0, 0, 0, 0}
			};
			c = getopt_long(argc, argv, "", long_options, &option_index);
			if ( c == -1)
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
						case 4: n_runs = atoi(optarg); break;
						case 5: n_episodes = atoi(optarg); break;
						case 6: n_steps = atoi(optarg); break;
						case 7: epsilon = atof(optarg); break; 
						case 8: N0 = atof(optarg); break;
						case 9: a = atof(optarg); break;
						case 10: environment_name = optarg; break;
						case 11: b_functions = atoi(optarg); break;
						case 12: grids = atoi(optarg); break;
						case 13: randomness = atof(optarg); break;
						case 14: episode_steps = atoi(optarg); break;
						case 15: n_samples = atoi(optarg); break;
						case 16: sampling = atoi(optarg); break;
						case 17: threshold = atof(optarg); break;
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
    assert (randomness >= 0 && randomness <= 1);
    assert (n_runs > 0);
    assert (n_episodes > 0);
    assert (n_steps > 0);
	assert (b_functions == 0 || b_functions == 1);
    assert (grid_size > 0);

    RandomNumberGenerator* rng;
	RandomNumberGenerator* environment_rng;

	MersenneTwisterRNG mersenne_twister_env;
    environment_rng = (RandomNumberGenerator*) &mersenne_twister_env;
	
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;
	
    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);
    environment_rng->manualSeed(228240153);
    rng->manualSeed(1361690241);
	
    std::cout << "Starting test program" << std::endl;
    
	std::cout << "Starting evaluation" << std::endl;
    // remember to use n_runs
    Statistics statistics;
    statistics.ep_stats.resize(n_episodes);
    statistics.reward.resize(n_steps,0.0);
    statistics.n_runs.resize(n_steps,0);
	
	
//	for (uint i=0; i<n_steps; ++i) {
//        statistics.reward[i] = 0;
//        statistics.n_runs[i] = 0;
//    }
	
	std::cout << " - Creating environment.." << std::endl;
	
	ContinuousStateEnvironment* environment = NULL;
	
	if (!strcmp(environment_name, "MountainCar")) {
		environment = new MountainCar();
		environment->setRandomness(randomness);
	}
	else if (!strcmp(environment_name, "Pendulum")) {
		environment = new Pendulum();
		environment->setRandomness(randomness);
	}
	else if (!strcmp(environment_name, "PuddleWorld")) {
		environment = new PuddleWorld();
	} else {
		fprintf(stderr, "Unknown environment %s \n", environment_name);
	}
	
	std::cout << environment_name << " environment creation completed" << std::endl;
	// making sure the number of states & actions is correct
	n_states	= environment->getNStates();
	n_actions	= environment->getNActions();
	
	std::cout <<  "Creating environment: " << environment_name
	<< " with " << n_states << " states and , "
	<< n_actions << " actions.\n";
	
	int m	= n_states + 1; //Input dimensions (input state dimension plus a dummy state)
	int d_r = 1;			// Reward dimensions
	int d_s = n_states;		//Output states dimensions 	
	
	RBFBasisSet* RBFs = NULL;
	if( b_functions == 1) {
		//Lower and upper environment bounds
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
	
		std::cout << "Creating Radial basis functions..." << std::endl;
		EvenGrid Discretisation(S_L, S_U, grids);
		RBFBasisSet* RBFs = new RBFBasisSet(Discretisation);
		
		m = RBFs->size() + 1; // redefinition of the input dimensions (size of the created basis functions plus a dummy state)

		std::cout << "# Number of basis functions: " << (RBFs->size() + 1)<< std::endl;
	}
	
	std::cout << "# Creating " << n_actions << " bayesian multivariate models..." << std::endl;
	
	//bayesian multivariate regression model for the system transition model
	std::vector<BayesianMultivariateRegression*> regression_t;
	regression_t.resize(n_actions);
	
	Matrix S0 = N0 * Matrix::Unity(d_s, d_s);
	for ( int i = 0; i < n_actions; ++i) {
		regression_t[i] = new BayesianMultivariateRegression(m, d_s, S0, N0, a);
	}
	
	//bayesian multivariate regression model for the system reward model
	std::vector<BayesianMultivariateRegression*> regression_r;
	regression_r.resize(n_actions);
	
	S0 = N0 * Matrix::Unity(d_r, d_r);
	for ( int i = 0; i < n_actions; ++i) {
		regression_r[i] = new BayesianMultivariateRegression(m, d_r, S0, N0, a);
	}
	
	std::cout << "Linear model initialization..." << std::endl;

	LinearModel<Vector,int>* lm = new LinearModel<Vector,int>( m, d_s, RBFs,environment);

	std::cout << "Representative model creation..." << std::endl;
	RepresentativeStateModel<LinearModel<Vector,int>, Vector, int> *RSM = new RepresentativeStateModel<LinearModel<Vector,int>,Vector,int>(gamma, threshold, *lm, n_samples, n_actions);
	
	
	std::cout << "Creating online algorithm" << std::endl;
	BayesianMultivariate* algorithm = NULL;
	algorithm = new BayesianMultivariate(n_actions,
										 m,
										 d_s,
										 gamma,
										 epsilon,
										 RBFs,
										 regression_t,
										 regression_r,
										 lm,
										 RSM);
	
	for ( uint run = 0; run < n_runs; ++run) {
		std::cout << "# Run: " << run << std::endl;

		TrainingAlgorithm(gamma, n_episodes, algorithm, environment, RSM);
		Statistics run_statistics = EvaluateAlgorithm(episode_steps,
                                                      n_episodes,
                                                      n_steps,
                                                      algorithm,
                                                      environment,
                                                      gamma);
		
		for (uint i=0; i<run_statistics.ep_stats.size(); ++i) {
			statistics.ep_stats[i].total_reward += run_statistics.ep_stats[i].total_reward;
			statistics.ep_stats[i].discounted_reward += run_statistics.ep_stats[i].discounted_reward;
			statistics.ep_stats[i].steps += run_statistics.ep_stats[i].steps;
			statistics.ep_stats[i].mse += run_statistics.ep_stats[i].mse;
			statistics.ep_stats[i].n_runs++;
		}
		
		for (uint i=0; i<run_statistics.reward.size(); ++i) {
			statistics.reward[i] += run_statistics.reward[i];
			statistics.n_runs[i]++;
		}
		algorithm->Reset();
	}

	delete environment;
	delete RBFs;
	for( int i = 0; i<n_actions; ++i){
		delete regression_t[i];
		delete regression_r[i];
	}
	delete lm;
	delete algorithm;
	
	for (uint i=0; i<statistics.ep_stats.size(); ++i) {
		statistics.ep_stats[i].total_reward /= (float) n_runs;
		statistics.ep_stats[i].discounted_reward /= (float) n_runs;
		statistics.ep_stats[i].steps /= n_runs;
		statistics.ep_stats[i].mse /= n_runs;
		std::cout << statistics.ep_stats[i].n_runs << " "
				  << statistics.ep_stats[i].total_reward << " "
				  << statistics.ep_stats[i].discounted_reward << " # EPISODE_RETURN"
				  << std::endl;
		std::cout << statistics.ep_stats[i].steps << " "
				  << statistics.ep_stats[i].mse << " # INST_PAYOFF"
				  << std::endl;
	}
	
	
	std::cout << "Done" << std::endl;
	
	return 0;
}

/*** Policy exploration
 n_steps: maximum number of total steps inside a trajectory.
 */

void TrainingAlgorithm(real gamma,
					   int n_episodes,
					   BayesianMultivariate* algorithm, 
					   ContinuousStateEnvironment* environment,
					   RepresentativeStateModel<LinearModel<Vector,int>, Vector, int> *R)
{
	std::cout << "#Training ...." << environment->Name() << std::endl;
	Vector state, next_state;
	environment->Reset();
	real reward;
	float epsilon_add = 1;
	int nsamples;
	bool action_ok = true;
	bool flag;
	
	for( int episode = 0; episode < n_episodes; ++episode)	{
		nsamples = R->getNSamples();
		printf("Number of representative states = %d\n",nsamples);
		for( int r_sample = 0; r_sample < nsamples; ++r_sample ) {
			environment->Reset();
			state = R->getSample(r_sample);
//			printf("Representative state (%d) = \n",r_sample); state.print(stdout);
			environment->setState(state);
			action_ok = true;
			int step = 0;
			flag = false;
			if(urandom() < epsilon_add) {
				flag = true;
			}
			while(action_ok && urandom() > (1 - gamma)) {
				step++;
				int action	= algorithm->Act(state);
				action_ok	= environment->Act(action);
				reward		= environment->getReward();
				next_state	= environment->getState();
				algorithm->Observe(state, action, reward, next_state);
				state = next_state;
				if(flag == true && urandom() < 0.1) {
					flag = false;
					R->AddState(next_state);
				}
			}
			//if(urandom() < epsilon_add) {
//				flag = true;
//			//	printf("ADD REPRESENTATIVE STATE\n");
////				R->AddState(next_state);
//			}
//			printf("Terminate at time step = %d\n",step);
		}
		epsilon_add = epsilon_add*0.5;
	
		algorithm->Update();
		printf ("# episode %d complete\n", episode);
	}
}
					   


/*** Evaluate an algorithm
 
 episode_steps: maximum number of steps per episode. If negative, then ignore
 n_steps: maximun number of total steps. If negative, then ignore.
 n_episodes: maximum number of episodes. Cannot be negative.
 */

Statistics EvaluateAlgorithm (int episode_steps,
                              int n_episodes,
                              uint n_steps,
                              BayesianMultivariate* algorithm,
                              ContinuousStateEnvironment* environment,
                              real gamma)
{
    std:: cout << "# evaluating..." << environment->Name() << std::endl;
    
	
    Statistics statistics;
    statistics.ep_stats.reserve(n_episodes); 
    statistics.reward.reserve(n_steps);
	
    real discount = 1.0;
    int current_time = 0;
    environment->Reset();

    std:: cout << "(running)" << std::endl;
    int episode = -1;
    bool action_ok = false;
    real total_reward = 0.0;
    real discounted_reward = 0.0;
	n_episodes = 1;
    for (uint step = 0; step < n_steps; ++step) {
        if (episode_steps > 0 && current_time >= episode_steps) {
            action_ok = false;
            environment->Reset();
        }
        if (!action_ok) {
            Vector state = environment->getState();
            real reward = environment->getReward();
            
			algorithm->Act(state);
		
            statistics.reward.resize(step + 1);
            statistics.reward[step] = reward;
            if (episode >= 0) {
                statistics.ep_stats[episode].steps++;
                statistics.ep_stats[episode].total_reward += reward;
				statistics.ep_stats[episode].discounted_reward += discount * reward;
            }
            total_reward += reward;
            discounted_reward += discount * reward;
			
            discount *= gamma;            
			
            episode++;
            statistics.ep_stats.resize(episode + 1);
            statistics.ep_stats[episode].total_reward = 0.0;
            statistics.ep_stats[episode].discounted_reward = 0.0;
            statistics.ep_stats[episode].steps = 0;
            discount = 1.0;
            environment->Reset();
          
            action_ok = true;
            current_time = 0;
            printf ("# episode %d complete\n", episode);
            if (n_episodes >= 0 && episode >= n_episodes) {
                logmsg (" Breaking after %d episodes,  %d steps\n",
						episode, step);
                break;
            }
            step++;
        }
		
        Vector state = environment->getState();
        real reward = environment->getReward();
        statistics.reward.resize(step + 1);
        statistics.reward[step] = reward;
        statistics.ep_stats[episode].steps++;
        statistics.ep_stats[episode].total_reward += reward;
        statistics.ep_stats[episode].discounted_reward += discount * reward;
        total_reward += reward;
        discounted_reward += discount * reward;
		
        discount *= gamma;
		
        int action;
		
		action = algorithm->Act(state);
           
        action_ok = environment->Act(action);
        current_time++;
		
    }
    printf(" %f %f # RUN_REWARD\n", total_reward, discounted_reward);
	
    if ((int) statistics.ep_stats.size() != n_episodes) {
        statistics.ep_stats.resize(statistics.ep_stats.size() - 1);
    }
    printf ("# Exiting after %d episodes, %d steps (%d %d)\n",
            episode, n_steps,
            (int) statistics.ep_stats.size(),
            (int) statistics.reward.size());
		
    return statistics;
}


#endif


