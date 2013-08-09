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
#include "CartPole.h"
#include "Pendulum.h"
#include "Acrobot.h"
#include "PuddleWorld.h"
#include "LinearDynamicQuadratic.h"
#include "DiscretisedEnvironment.h"

// -- Randomness -- //
#include "RandomNumberFile.h"
#include "MersenneTwister.h"

// -- Algorithms and models --//
#include "Grid.h"
#include "BasisSet.h"
#include "CoverTree.h"
#include "CoverBayesianMultivariate.h"
#include "CoverFittedValueIteration.h"
#include "CoverFittedLSTD.h"
#include "CoverFittedLSTDQ.h"


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

Statistics TrainingAlgorithm(real gamma,
							  int episode_steps,
							  int n_episodes,
							  CoverBayesianMultivariate* algorithm, 
							  ContinuousStateEnvironment* environment);

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
--environment:		{PuddleWorld, MountainCar, Pendulum, Linear, Acrobot, CartPole} \n\
--gamma:			reward discounting in [0,1] (* 0.99)\n\
--randomness:		environment randomness (* 0.0)\n\
--n_runs:			maximum number of runs (* 1)\n\
--episodes:			maximum number of train episodes (ignored if < 0)\n\
--episode_steps:    maximum number of steps in each episode (ignored if <0)\n\
--rbf:				use of basis functions or not (* 0 or 1)\n\
--grids:			number of grid intervals for discretised environments (*3)\n\
--epsilon:			use epsilon-greedy with randomness in [0,1] (* 0.01)\n\
--n_samples:		number of samples that used on the fitting (*3000)\n\
--a:				first linear model parameter (*0.1) \n\
--N0:				second linear model parameter (*0.1) \n\
--scale: rbf scale (1)\n\
\n\
* denotes default parameter\n\
\n";


int main(int argc, char* argv[])
{
	real gamma					= 0.999;
	real randomness				= 0.0;
	real epsilon				= 0.0;
	real tree_c					= 1.0;
	real N0						= 1.0;
	real a						= 1.0;
	uint rbf					= 0;
	uint grids					= 4;
	real scale					= 1.0;
	int episodes				= 1000;
	int episode_steps			= 1000;
	int n_samples				= 3000;
	uint n_runs					= 1;
	
	ulong seed = time(NULL);
	char* seed_filename			= 0;

	
	
	const char * environment_name = "MountainCar";
	{
		//options
		int c;
		int digit_optind = 0;
		while(1) {
			int this_option_optind = optind ? optind : 1;
			int option_index = 0;
			static struct option long_options[] = {
				{"gamma", required_argument, 0, 0},						//0
				{"randomness", required_argument, 0, 0},				//1
                {"epsilon", required_argument, 0, 0},					//2
				{"tree_c", required_argument, 0, 0},					//3
				{"N0", required_argument, 0, 0},						//4
				{"a", required_argument, 0, 0},							//5
				{"rbf", required_argument, 0, 0},						//6
                {"grids", required_argument, 0, 0},						//7
				{"scale", required_argument, 0, 0},						//8
				{"episodes", required_argument, 0, 0},					//9
				{"episode_steps", required_argument, 0, 0},				//10
				{"environment_name", required_argument, 0, 0},			//11
				{"n_samples", required_argument, 0, 0},					//12
				{"n_runs", required_argument, 0, 0},					//13
				{"seed_filename", required_argument, 0, 0},				//14
				{0, 0, 0, 0}
			};
			c = getopt_long(argc, argv, "", long_options, &option_index);
			if ( c == -1)
				break;
			
			switch (c) {
				case 0:
#if 1
					printf ("option %s (%d)", long_options[option_index].name, option_index);
					if (optarg)
						printf (" with arg %s", optarg);
					printf ("\n");
#endif
					switch (option_index) {
						case 0: gamma = atof(optarg); break;
						case 1: randomness = atof(optarg); break;
						case 2: epsilon = atof(optarg); break; 
						case 3: tree_c = atof(optarg); break;
						case 4: N0 = atof(optarg); break;
						case 5: a = atof(optarg); break;
						case 6: rbf = atoi(optarg); break;
						case 7: grids = atoi(optarg); break;
						case 8: scale = atof(optarg); break;
						case 9: episodes = atoi(optarg); break;
						case 10: episode_steps = atoi(optarg); break;
						case 11: environment_name = optarg; break;
						case 12: n_samples = atoi(optarg); break;
						case 13: n_runs = atoi(optarg); break;
						case 14: seed_filename = optarg; break;
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
	
	assert(gamma >= 0 && gamma <=1);
	assert(randomness >=0 && randomness <= 1);
	assert(epsilon >= 0);
	assert(N0 > 0);
	assert(a > 0);
	assert(rbf == 0 || rbf == 1);
	assert(grids > 1);
	assert(scale > 0);
	assert(episodes > 0);
	assert(episode_steps > 0);
	assert(n_samples > 0);
		
//  RandomNumberGenerator* rng;
//	RandomNumberGenerator* environment_rng;
//	
//	MersenneTwisterRNG mersenne_twister_env;
//  environment_rng = (RandomNumberGenerator*) &mersenne_twister_env;
//	
//  MersenneTwisterRNG mersenne_twister;
//  rng = (RandomNumberGenerator*) &mersenne_twister;
//	
//  srand48(34987235);
//  srand(34987235);
//  setRandomSeed(34987235);
//  environment_rng->manualSeed(228240153);
//  rng->manualSeed(1361690241);
	
	MersenneTwisterRNG rng;
   	
	if (seed_filename) {
        RandomNumberFile rnf(seed_filename);
        rnf.manualSeed(seed);
        seed = rnf.random();
    }
	
    logmsg("seed: %ld\n", seed);
    srand(seed);
    srand48(seed);
    rng.manualSeed(seed);
    setRandomSeed(seed);
	
    std::cout << "Starting test program" << std::endl;
    
	std::cout << "Starting evaluation" << std::endl;
	
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
	else if (!strcmp(environment_name, "CartPole")) {
		environment = new CartPole();
		environment->setRandomness(randomness);
	}
	else if (!strcmp(environment_name, "Acrobot")) {
		environment = new Acrobot();
		environment->setRandomness(randomness);
	}
	else if (!strcmp(environment_name, "PuddleWorld")) {
		environment = new PuddleWorld();
	}
	else if(!strcmp(environment_name,"Linear")) {
		environment = new LinearDynamicQuadratic();
	}
	else {
		fprintf(stderr, "Unknown environment %s \n", environment_name);
	}
	
	std::cout << environment_name << " environment creation completed" << std::endl;
	// making sure the number of states & actions is correct
	int n_states	= environment->getNStates();
	int n_actions	= environment->getNActions();
	
	std::cout <<  "Creating environment: " << environment_name
	<< " with " << n_states << " states and , "
	<< n_actions << " actions.\n";
	
	int m	= n_states + 1; //Input dimensions (input state dimension plus a dummy state)
	
	RBFBasisSet* RBFs = NULL;
	
	if( rbf == 1) {
		Vector S_L	= environment->StateLowerBound();
		S_L.print(stdout);
		Vector S_U	= environment->StateUpperBound();
		S_U.print(stdout);		
		//Lower and upper environment bounds
		std::cout << "Creating Radial basis functions..." << std::endl;
		EvenGrid Discretisation(S_L, S_U, grids);
		
		RBFs = new RBFBasisSet(Discretisation,1.0);
		
		m = RBFs->size() + 1; // redefinition of the input dimensions (size of the created basis functions plus a dummy state)
		std::cout << "# Number of basis functions: " << m << std::endl;
	}
	
	
	//Cover Trees creation (one for each action).
	std::cout << "--Creating environment models (one CoverTree for each action)..." << std::endl;
	std::vector<CoverTree*> cover;
	cover.resize(n_actions);
	
	for( int i = 0; i < n_actions; ++i) {
		cover[i] = new CoverTree(tree_c, a, N0, RBFs);
	}
	std::cout << n_actions << " Cover trees was created..." << std::endl;
	
	std::cout << "Cover Fitted value iteration initialization..." << std::endl;
	CoverFittedValueIteration<Vector,int> *FVI = NULL;
	//	CoverFittedValueIteration<Vector,int> *FVI = new CoverFittedValueIteration<Vector,int>(gamma, 2000, 5, environment, cover);
	
	std::cout << "Fitted LSTD initialization..." << std::endl;
	CoverFittedLSTD<Vector,int> *FLSTD = new CoverFittedLSTD<Vector,int>(gamma, n_samples, 1, environment, cover, grids, scale);
	//	CoverFittedLSTD<Vector,int> *FLSTD = new CoverFittedLSTD<Vector,int>(gamma, 3000, 1, environment, cover);
	//	CoverFittedLSTD<Vector,int> *FLSTD = NULL;
	
	//	CoverFittedLSTDQ<Vector,int> *FLSTDQ = new CoverFittedLSTDQ<Vector,int>(gamma, 5000, 1, environment, cover, 4, 1.0);
	//	CoverFittedLSTDQ<Vector,int>* FLSTDQ = new CoverFittedLSTDQ<Vector,int>(gamma, 3000, 1, environment, cover);
	CoverFittedLSTDQ<Vector,int> *FLSTDQ = NULL;

	Matrix Stats(n_runs, episodes);
	Matrix Statr(n_runs, episodes);
	
	std::cout << "Basic algorithm creation" << std::endl;
	CoverBayesianMultivariate* algorithm = new CoverBayesianMultivariate( n_actions,
																		 gamma,
																		 epsilon,
																		 cover,
																		 FVI,
																		 FLSTD,
																		 FLSTDQ);
	
//	for ( uint run = 0; run < n_runs; ++run) {
	uint run = 0;
	std::cout << "# Run: " << run << std::endl;
	//		algorithm->setGeometricSchedule(0.01,0.1);
	//		TrainingAlgorithm1(gamma, n_train_episodes, algorithm, environment, RSM);
		
	Statistics run_statistics = TrainingAlgorithm(gamma, episode_steps, episodes, algorithm, environment);
		
	epsilon = 0.0;
		
	real train_steps = 0.0;
		
	for( int i = 0; i < episodes; ++i) {
		Stats(run,i) = run_statistics.ep_stats[i].steps;
		Statr(run,i) = run_statistics.ep_stats[i].total_reward;
		train_steps += run_statistics.ep_stats[i].steps;
	}
	printf("Mean number of steps = %f\n", (train_steps / episodes));
		
	algorithm->Reset();
//	}
	
	char buffer[100];
	sprintf (buffer, "CLBRL_ONLINE_RESULTS_STEPS_%s",environment_name);
	FILE *output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Stats.print(output);
	}
	fclose(output);
	sprintf (buffer, "CLBRL_ONLINE_RESULTS_REWARDS_%s",environment_name);
	output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Statr.print(output);
	}
	fclose(output);
	
	//Pointer clearness
	delete environment;
	delete RBFs;
	
	for( int i = 0; i<n_actions; ++i){
		delete cover[i];
	}
	delete algorithm;
	
	std::cout << "Done" << std::endl;
	
	return 0;
}

/*** Policy exploration
 n_steps: maximum number of total steps inside a trajectory.
 */
Statistics TrainingAlgorithm(real gamma,
							  int episode_steps,
							  int n_episodes,
							  CoverBayesianMultivariate* algorithm, 
							  ContinuousStateEnvironment* environment)
{
	std::cout << "#Training ...." << environment->Name() << std::endl;
	Vector state;
	Vector next_state;
	environment->Reset();
	real reward;
	real discount = gamma;
	int n_actions = environment->getNActions(); //Possible actions.
	
	///Statistics.
    Statistics statistics;
	statistics.ep_stats.reserve(n_episodes); 
	int step = 0;
	int initial_iter = 2;
	//	int nsamples;
	bool action_ok = true;
	int action;
	
	statistics.ep_stats[0].total_reward = 0.0;
	statistics.ep_stats[0].discounted_reward = 0.0;
	statistics.ep_stats[0].steps = 0;
	
	for( int episode = 0; episode < n_episodes; ++episode)	{
		
		environment->Reset();
		state   = environment->getState();
//		printf("Initial State => ");
//		state.print(stdout);
		if(episode<initial_iter) {
			action = (int) floor(urandom(0.0, (real) n_actions));  //Randomly selected action
		}
		else {
			action	= algorithm->Act(state);
		}

		action_ok = true;
		step = 0;
		
		while(action_ok && step < episode_steps) {
			step++;
			action_ok	= environment->Act(action);
			
			reward		= environment->getReward();
			next_state	= environment->getState();
			
			algorithm->Observe(state, action, reward, next_state);
			
			state = next_state;
			statistics.ep_stats[episode].steps++;
			statistics.ep_stats[episode].total_reward += reward;
			statistics.ep_stats[episode].discounted_reward += discount * reward;
			discount *= gamma;
//			action	= algorithm->Act(state);
			if(episode < initial_iter) {
				action = (int) floor(urandom(0.0, (real) n_actions));  //Randomly selected action
			}
			else {
				action	= algorithm->Act(state);
			}
		}
		printf ("# episode %d complete -> Steps = %d\n", episode, step);
		if(episode>=(initial_iter-1)) {
			algorithm->Update();
		}
		///Statistics Initializastion
		if((episode + 1) < n_episodes) {
			statistics.ep_stats[episode + 1].total_reward = 0.0;
			statistics.ep_stats[episode + 1].discounted_reward = 0.0;
			statistics.ep_stats[episode + 1].steps = 0;
		}
	}
	
	//algorithm->Update();
	//algorithm->Predict(states);
	return statistics;
}

#endif


