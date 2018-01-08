/* -*- Mode: C++; -*- */
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
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
#include "Acrobot.h"
#include "Gridworld.h"
#include "OneDMaze.h"
#include "DiscreteChain.h"
#include "RiverSwim.h"
#include "CartPole.h"
#include "Bike.h"
#include "OptimisticTask.h"
#include "InventoryManagement.h"
#include "DoubleLoop.h"

// -- Continuous environments -- //
#include "MountainCar.h"
#include "Pendulum.h"
#include "PuddleWorld.h"
#include "LinearDynamicQuadratic.h"
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
#include "FittedQValueIteration.h"
#include "FittedValueIteration.h"
#include "FittedLSTD.h"
#include "BayesianMultivariate.h"


// -- Usual includes -- //
#include <cstring>
#include <getopt.h>
/** Options */
struct Options
{
    real gamma;						///< discount factor
    const char* environment_name;	///< environment name
	const char* method;				///< Method API
	const char* sampling;			///< Sampling Approach (Marginal, Thompson)
	real a_model;					///< Bayesian Linear Regression model parameter a
	real N0_model;					///< Bayesian Linear Regression model parameter N0
    RandomNumberGenerator& rng;		///< random number generator
    int n_samples;					///< number of sampled pairs for LSTD/FVI training
	int n_training_e_steps;			///< Max length of training trajectories
    int n_training_episodes;		///< number of training trajectories
	int n_testing_e_steps;			///< Max length of testing trajectories.
    int n_testing_episodes;			///< number of testing trajectories
    int grid;						///< grid size for features
    real scale;						///< scale of RBF basis
	real randomness;				///< environment randomness
	real epsilon;					///< E-greedy action selection probability
	int n_runs;					    ///< How many times we will run the experiment
	
    Options(RandomNumberGenerator& rng_) :
	gamma(0.999),
	environment_name("MountainCar"),
	method("LSTD"),
	sampling("Thompson"),
	a_model(0.001),
	N0_model(0.001),
	rng(rng_), 
	n_samples(5000),
	n_training_e_steps(40),
	n_training_episodes(1000),
	n_testing_e_steps(1000),
	n_testing_episodes(1000),
	grid(4),
	scale(1.0),
	randomness(0.0),
	epsilon(0.0),
	n_runs(100)
    {
    }
	
    void ShowOptions()
    {
        logmsg("---------------------------\n");
		logmsg("===========================\n");
        logmsg("Options\n");
        logmsg("===========================\n");
        logmsg("Discount: %f\n", gamma);
		logmsg("Environment: %s\n", environment_name);
		logmsg("Method: %s\n", method);
		logmsg("Sampling: %s\n", sampling);
		logmsg("a parameter: %f\n", a_model);
		logmsg("N0 parameter: %f\n", N0_model);
        logmsg("n_samples: %d\n", n_samples);
        logmsg("n_training_e_steps: %d\n",  n_training_e_steps);
		logmsg("n_training_episodes: %d\n", n_training_episodes);
		logmsg("n_testing_e_steps: %d\n", n_testing_e_steps);
		logmsg("n_testing_episodes: %d\n", n_testing_episodes);
        logmsg("grid %d\n", grid);
        logmsg("scale %f\n", scale);
		logmsg("epsilon %f\n", epsilon);
		logmsg("num of runs %d", n_runs);
        logmsg("---------------------------\n");
    }
};

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

void OfflineLBRL(real gamma,
				 int n_steps, 
				 int n_episodes,
				 BayesianMultivariate* algorithm, 
				 ContinuousStateEnvironment* environment);

Statistics OnlineLBRL(real gamma,
					  int n_steps, 
					  int n_episodes,
					  BayesianMultivariate* algorithm, 
					  ContinuousStateEnvironment* environment);

void EvaluatePrediction(int N,
						BayesianMultivariate* algorithm,
						ContinuousStateEnvironment* environment);

Statistics EvaluateAlgorithm(real gamma,
							int episode_steps,
							int n_episodes,
							BayesianMultivariate* algorithm,
							ContinuousStateEnvironment* environment);

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
--environment:			{MountainCar, Pendulum, Puddle, Bicycle, CartPole, Acrobot}\n\
--method:				{FVI,LSTD,FQVI,RSM}\n\
--sampling:				{Thompson, Marginal} (*Thompson)\n\
--a_model:				enviroment model parameter a \n\
--N0_model:				enviroment model parameter N0 \n\
--discount:				reward discounting in [0,1]\n\
--n_samples:			number of selected (uniformly random) samples that used in the policy optimization scheme (LSTD)\n\
--n_training_e_steps:	maximum horizon of training rollouts (episodes)\n\
--n_training_episodes:  number of training episodes\n\
--n_testing_e_steps:	maximum horizon of testing episodes\n\
--n_testing_episodes:	number of testing episodes\n\
--grid:                 number of grid intervals for LSTD\n\
--scale:                RBF scale for LSTD\n\
--online:               do the online test\n\
--seed:                 seed all the RNGs with this\n\
--seed_file:            select a binary file to choose seeds from (use in conjunction with --seed to select the n-th seed in the file)\n\
--randomness:			environment randomness (* 0.0)\n\
--epsilon:				e-greedy action selection probability (*0.0)\n\
--n_runs:				the number of times where we will run the experiment (*100)\n\
\n\
* denotes default parameter\n\
\n";


int main(int argc, char* argv[])
{
	ulong seed = time(NULL);
	char* seed_filename = 0;
	const char* test_mode = "offline";
	bool online_test = false;
	bool model_rbf = false;
	
	MersenneTwisterRNG rng;
	Options options(rng);
	{
		//options
		int c;
		int digit_optind = 0;
		while(1) {
			int this_option_optind = optind ? optind : 1;
			int option_index = 0;
			static struct option long_options[] = {
				{"discount", required_argument, 0, 0},				//0
				{"environment", required_argument, 0, 0},			//1
				{"method", required_argument, 0, 0},				//2
				{"sampling", required_argument, 0, 0},				//3
				{"a_model", required_argument, 0, 0},				//4
				{"N0_model", required_argument, 0, 0},				//5
				{"n_samples", required_argument, 0, 0},				//6
				{"n_training_e_steps", required_argument, 0, 0},	//7
				{"n_training_episodes", required_argument, 0, 0},	//8
				{"n_testing_e_steps", required_argument, 0, 0},		//9
				{"n_testing_episodes", required_argument, 0, 0},	//10
				{"grid", required_argument, 0, 0},					//11
				{"scale", required_argument, 0, 0},					//12
				{"online", no_argument, 0, 0},						//13
				{"seed", required_argument, 0, 0},					//14
				{"seed_file", required_argument, 0, 0},				//15
				{"randomness", required_argument, 0, 0},			//16
				{"epsilon", required_argument, 0, 0},				//17
				{"rbf_model", required_argument, 0, 0},				//18
				{"n_runs", required_argument, 0, 0},				//19
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
						case 0: options.gamma = atof(optarg); break;
						case 1: options.environment_name = optarg; break;
						case 2: options.method = optarg; break;
						case 3: options.sampling = optarg; break;
						case 4: options.a_model = atof(optarg); break;
						case 5: options.N0_model = atof(optarg); break;
						case 6: options.n_samples = atoi(optarg); break;
						case 7: options.n_training_e_steps = atoi(optarg); break;
						case 8: options.n_training_episodes = atoi(optarg); break;
						case 9: options.n_testing_e_steps = atoi(optarg); break;
						case 10: options.n_testing_episodes = atoi(optarg); break;
						case 11: options.grid = atoi(optarg); break;
						case 12: options.scale = atof(optarg); break;
						case 13: online_test = true; break;
						case 14: seed = atoi(optarg); break;
						case 15: seed_filename = optarg; break;
						case 16: options.randomness = atof(optarg); break;
						case 17: options.epsilon = atof(optarg); break;
						case 18: model_rbf = true; break;
						case 19: options.n_runs = atoi(optarg); break;
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
	
	if (options.gamma < 0 || options.gamma > 1) {
        Serror("Discount factor gamma must be in [0,1]\n");
        exit(-1);
    }
	
	if (!options.environment_name) {
        Serror("Must specify environment\n");
        exit(-1);
    }
	
	if (!options.method) {
        Serror("Must specify API method\n");
        exit(-1);
    }
	
	if (options.a_model < 0) {
		Serror("\alpha parameter must be > 0\n");
        exit(-1);
    }
	
	if (options.N0_model < 0) {
		Serror("N0 parameters must be > 0\n");
        exit(-1);
    }	
	
    if (options.n_samples < 1) {
        Serror("n_samples must be >= 1\n");
        exit(-1);
    }
	
	if (options.n_training_e_steps < 1) {
        Serror("n_training_e_steps must be >= 1\n");
        exit(-1);
    } 
	
    if (options.n_training_episodes < 1) {
        Serror("n_training_episodes must be >= 1\n");
        exit(-1);
    }
	
	if (options.n_testing_e_steps < 1) {
        Serror("n_training_e_steps must be >= 1\n");
        exit(-1);
    } 
    
	if (options.n_testing_episodes < 0) {
        Serror("n_testing_episodes must be >= 1\n");
        exit(-1);
    }
	
	if(options.scale < 0) {
		Serror("scale must be > 0\n");
		exit(-1);
	}	
	if(options.n_runs < 0) {
		Serror("the number of runs must be positive\n");
		exit(-1);
	}
	
   	if(seed_filename) {
		RandomNumberFile rnf(seed_filename);
		rnf.manualSeed(seed);
		seed = rnf.random();
	}
	
	logmsg("seed: %ld\n", seed);
	srand(seed);
	srand48(seed);
	rng.manualSeed(seed);
	setRandomSeed(seed);
	
    std::cout << "Starting GPRL" << std::endl;
   	
	options.ShowOptions();
	
	std::cout << " - Creating environment.." << std::endl;
	
	ContinuousStateEnvironment* environment = NULL;
	
	if (!strcmp(options.environment_name, "MountainCar")) {
		environment = new MountainCar();
	}
	else if (!strcmp(options.environment_name, "Pendulum")) {
		environment = new Pendulum();
	}
	else if (!strcmp(options.environment_name, "PuddleWorld")) {
		environment = new PuddleWorld();
	}
	else if (!strcmp(options.environment_name, "Acrobot")) {
		environment = new Acrobot();
	}
	else if (!strcmp(options.environment_name, "CartPole")) {
		environment = new CartPole();
	}
	else if (!strcmp(options.environment_name, "Bike")) {
		environment = new Bike();
	}
	else if(!strcmp(options.environment_name,"Linear")) {
		environment = new LinearDynamicQuadratic();
	}
	else {
		fprintf(stderr, "Unknown environment %s \n", options.environment_name);
	}
	environment->setRandomness(options.randomness);
	
	int n_states	= environment->getNStates();      /// Environment dimensions
	int n_actions	= environment->getNActions();     /// Number of environment's actions
	Vector S_L		= environment->StateLowerBound(); /// Environment lower bounds
	Vector S_U		= environment->StateUpperBound(); /// Environment upper bounds
	
	std::cout <<  "Creating environment: " << options.environment_name
	<< " with " << n_states << " states and , "
	<< n_actions << " actions.\n";
	
	int m	= n_states + 1; //Input dimensions (input state dimension plus a dummy state)
	int d_r = 1;			//Reward dimensions
	int d_s = n_states;		//Output states dimensions 	
	srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);
	RBFBasisSet* RBFs = NULL;
	
	if( model_rbf ) {
		Vector S_L	= environment->StateLowerBound();
		S_L.print(stdout);
		Vector S_U	= environment->StateUpperBound();
		S_U.print(stdout);		
		//Lower and upper environment bounds
		std::cout << "Creating Radial basis functions..." << std::endl;
		EvenGrid Discretisation(S_L, S_U, options.grid);
		
		RBFs = new RBFBasisSet(Discretisation,1.0);
		
		m = RBFs->size() + 1; //redefinition of the input dimensions (size of the created basis functions plus a dummy state)
		std::cout << "# Number of basis functions: " << m << std::endl;
	}
	
	std::cout << "# Creating " << n_actions << " bayesian multivariate models..." << std::endl;
	
	//Bayesian Multivariate Regression model for the system transition model, one for each action
	std::vector<BayesianMultivariateRegression*> regression_t(n_actions);
	//Bayesian Multivariate Regression model for the system reward model, one for each action
	std::vector<BayesianMultivariateRegression*> regression_r(n_actions);

	Matrix S0s = options.N0_model * Matrix::Unity(d_s, d_s);
	Matrix S0r = options.N0_model * Matrix::Unity(d_r, d_r);
	for ( int i = 0; i < n_actions; ++i) {
		if(!strcmp(options.sampling, "Thompson")) {
			regression_t[i] = new BayesianMultivariateRegression(m, d_s, S0s, options.N0_model, options.a_model, true);
			regression_r[i] = new BayesianMultivariateRegression(m, d_r, S0r, options.N0_model, options.a_model, true);
		} else if(!strcmp(options.sampling, "Marginal")) {
			regression_t[i] = new BayesianMultivariateRegression(m, d_s, S0s, options.N0_model, options.a_model, false);
			regression_r[i] = new BayesianMultivariateRegression(m, d_r, S0r, options.N0_model, options.a_model, false);		
		} else {
			fprintf(stderr, "Unknown Sampling approach %s\n", options.sampling);
		}
		
	}
	
	std::cout << "Creation of the bayesian multivariate models completed..." << std::endl;
	
	std::cout << "API Method initialization" << std::endl;
	LinearModel<Vector,int>* lm = NULL;
	RepresentativeStateModel<LinearModel<Vector,int>, Vector, int> *RSM = NULL;
	FittedValueIteration<Vector,int> *FVI = NULL;
	FittedValueIteration<Vector,int> *private_FVI = NULL;
	FittedQValueIteration<Vector,int> *FQVI = NULL;
	FittedLSTD<Vector,int> *FLSTD = NULL;
	
	if (!strcmp(options.method, "FVI")) {
		std::cout << "Fitted Value Iteration (FVI) Initialization" << std::endl;
		FVI = new FittedValueIteration<Vector,int>(options.gamma, options.n_samples, 1, options.grid, environment, regression_t, RBFs, options.scale);
	}
	else if(!strcmp(options.method,"LSTD")) {
		std::cout << "Least Square Temporal Difference (LSTD) Initialization" << std::endl;
		FLSTD = new FittedLSTD<Vector,int>(options.gamma, options.n_samples, 1, options.grid, environment, regression_t, RBFs, options.scale);
	} else if (!strcmp(options.method, "FQVI")) {
		std::cout << "Fitted Q-Value Iteration (FQVI) Initialization" << std::endl;
		FQVI = new FittedQValueIteration<Vector,int>(options.gamma, options.n_samples, 1, options.grid, environment, regression_t, RBFs, options.scale);
	} else if (!strcmp(options.method, "RSM")) {
		std::cout << "Representative State Model (RSM) Initialization" << std::endl;
		private_FVI = new FittedValueIteration<Vector,int>(options.gamma, options.n_samples, 1, options.grid, environment, regression_t, RBFs, options.scale);
		RSM = NULL;
		//new RepresentativeStateModel<LinearModel<Vector,int>,Vector,int>(gamma, 0.000001, *lm, n_states, n_actions, private_FVI);
	} else {
		fprintf(stderr, "Unknown method %s \n", options.method);
	}
	
	BayesianMultivariate* algorithm = new BayesianMultivariate(n_actions,
															   m,
															   d_s,
															   options.gamma,
															   options.epsilon,
															   RBFs,
															   regression_t,
															   regression_r,
															   lm,
															   RSM,
															   FVI,
															   FLSTD,
															   FQVI);
	
	Statistics run_statistics;
	Matrix Stats(options.n_runs,options.n_testing_e_steps);
	Matrix Statr(options.n_runs,options.n_testing_episodes);
	for(int n = 0; n<options.n_runs; ++n) {
		int episodes;
		if(online_test) {	
			test_mode ="online";
			episodes = options.n_training_episodes;
			run_statistics = OnlineLBRL(options.gamma,
										options.n_testing_e_steps, 
										options.n_testing_episodes,
										algorithm,
										environment);
		} else {
			episodes = options.n_testing_episodes;
			OfflineLBRL(options.gamma,
						options.n_training_e_steps,
						options.n_training_episodes,
						algorithm,
						environment);
		
			run_statistics = EvaluateAlgorithm(options.gamma,
										   options.n_testing_e_steps,
										   options.n_testing_episodes,
										   algorithm, 
										   environment);
		}
	
	
		real train_steps = 0.0;
		for( uint i = 0; i < run_statistics.ep_stats.size(); ++i) {
			Stats(n, i) = run_statistics.ep_stats[i].steps;
			Statr(n, i) = run_statistics.ep_stats[i].total_reward;
			train_steps += run_statistics.ep_stats[i].steps;
		}
		printf("Mean number of steps = %f\n", train_steps / episodes);
		algorithm->Reset();
	}
	char buffer[100];
	sprintf (buffer, "%s_LBRL_%s_%d_%s_%s.steps",
			 test_mode,
			 options.sampling,
			 options.grid,
			 options.method,
			 options.environment_name);
	FILE *output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Stats.print(output);
	}
	fclose(output);
	
	sprintf (buffer, "%s_LBRL_%s_%d_%s_%s.reward",
			 test_mode,
			 options.sampling,
			 options.grid,
			 options.method,
			 options.environment_name);
	output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Statr.print(output);
	}
	fclose(output);
	
	for( int i = 0; i<n_actions; ++i){
			delete regression_t[i];
			delete regression_r[i];
	}
	delete RBFs;
	delete FVI;
	delete private_FVI;
	delete RSM;
	delete FLSTD;
	delete algorithm;
	delete environment;
	
	std::cout << "Done" << std::endl;
	
	return 0;
}

void OfflineLBRL(real gamma,
				 int n_steps, 
				 int n_episodes,
				 BayesianMultivariate* algorithm, 
				 ContinuousStateEnvironment* environment)
{
	std::cout << "#Offline Training ...." << environment->Name() << std::endl;
	Vector state;
	Vector next_state;
	real reward;
	bool action_ok;
	int current_time, action;

	for( int episode = 0; episode < n_episodes; ++episode)	{
		environment->Reset();
		state = environment->getState();
		action = (int) floor(urandom(0.0, (real) environment->getNActions()));
		
		current_time = 0;
		action_ok = true;
		
		while(action_ok && current_time < n_steps) {
			action_ok	= environment->Act(action);
			
			reward		= environment->getReward();
			next_state	= environment->getState();

			algorithm->Observe(state, action, reward, next_state);
			
			state = next_state;
			
			action  = (int)floor(urandom(0.0, (real) environment->getNActions()));
			current_time++;
		}
	}
	std::vector<Vector> states;

//	EvaluatePrediction(5000, algorithm, environment);
	
	algorithm->Update();
	
	printf("The offline LBRL training was completed\n");
}

Statistics OnlineLBRL(real gamma,
				int n_steps, 
				int n_episodes,
				BayesianMultivariate* algorithm, 
				ContinuousStateEnvironment* environment)
{
	std::cout << "#Online training..." << environment->Name() << std::endl;
	Vector state, next_state;
	int action;
	real reward;
	int n_actions = environment->getNActions();
	int n_states  = environment->getNStates();
	std::vector< std::vector<Vector> > Input(n_actions);
	std::vector< std::vector<real> > Output(n_actions*n_states);
	
	///Statistics.
    Statistics statistics;
	statistics.ep_stats.reserve(n_episodes); 
	statistics.reward.reserve(n_steps*n_episodes);
	
	statistics.ep_stats.resize(1);
	statistics.ep_stats[0].total_reward = 0.0;
	statistics.ep_stats[0].discounted_reward = 0.0;
	statistics.ep_stats[0].steps = 0;
	
	for(int episode = 0; episode <n_episodes; ++episode) {
		environment->Reset();
		state	= environment->getState();
		action	= algorithm->Act(state);
		
		real discount	= 1;
		int step		= 0;
		bool action_ok	= true;
		
		while(action_ok && step < n_steps) {
			step++;
			action_ok	= environment->Act(action);
			reward		= environment->getReward();
			
			next_state  = environment->getState();
			
			algorithm->Observe(state, action, reward, next_state);
			
			action = algorithm->Act(next_state);
			state = next_state;
			
			statistics.ep_stats[episode].steps++;
			statistics.ep_stats[episode].total_reward += reward;
			statistics.ep_stats[episode].discounted_reward += discount * reward;
			discount *= gamma;
		}
		///Statistics Initializastion
		if((episode+2) < n_episodes+1) {
			statistics.ep_stats.resize(episode+2);
			statistics.ep_stats[episode+1].total_reward = 0.0;
			statistics.ep_stats[episode+1].discounted_reward = 0.0;
			statistics.ep_stats[episode+1].steps = 0;
		}
		printf ("# episode %d complete -> Steps = %d\n", episode, step);
		algorithm->Update();
	}
	printf("The offline LBRL training was completed\n");
	
	return statistics;
}

void EvaluatePrediction(int N,
						BayesianMultivariate* algorithm,
						ContinuousStateEnvironment* environment)
{
	Vector S_L	= environment->StateLowerBound();
	Vector S_U	= environment->StateUpperBound();

	std::vector<Vector> states;
	for(int i=0; i < N; ++i) {
		states.push_back(urandom(S_L, S_U));
	}
	FILE *input = fopen("Input_samples.txt","w");
	if(input!=NULL) {
		for( uint i = 0; i < states.size(); ++i) {
			states[i].print(input);
		}
	}
	fclose(input);
	char buffer[100];
	int n;
	for ( int a = 0; a < (int)environment->getNActions(); ++a) {
		n = sprintf (buffer, "Output_samples_action_%d", a);
		FILE *output	= fopen(buffer,"w");
		n = sprintf (buffer, "Reward_samples_action_%d", a);
		FILE *rew		= fopen(buffer,"w");
		if(output!=NULL) {
			for( uint i = 0; i < states.size(); ++i) {
				environment->Reset();
//				state		= states[i];
				environment->setState(states[i]);
				environment->Act(a);
				real reward	= environment->getReward();
				Vector r(reward);
//				r.print(rew);
				Vector next_state	= environment->getState();
				next_state.print(output);
			}
		}
		fclose(output);
		fclose(rew);
	}
	algorithm->Predict(states);
}

/*** Evaluate an algorithm
 episode_steps: maximum number of steps per episode. If negative, then ignore
 n_steps: maximun number of total steps. If negative, then ignore.
 n_episodes: maximum number of episodes. Cannot be negative.
 */
Statistics EvaluateAlgorithm(real gamma,
							  int episode_steps,
                              int n_episodes,
                              BayesianMultivariate* algorithm,
                              ContinuousStateEnvironment* environment)
{
    std:: cout << "# evaluating... " << environment->Name() << std::endl;
	int n_steps = episode_steps*n_episodes;
    Statistics statistics;
    statistics.ep_stats.reserve(n_episodes); 
    statistics.reward.reserve(n_steps);
	uint step = 0;
	
    real discount = gamma;
    int current_time = 0;
    environment->Reset();
	
	//	environment->getState().print(stdout);
    std:: cout << "(running)" << std::endl;
    int episode = -1;
    bool action_ok = false;
    real total_reward = 0.0;
    real discounted_reward = 0.0;
    
	while(1) {
		step++;
        if (episode_steps > 0 && current_time >= episode_steps) {
            action_ok = false;
            environment->Reset();
        }
        if (!action_ok) {
            Vector state	= environment->getState();
            real reward		= environment->getReward();
			
            statistics.reward.resize(step + 1);
            statistics.reward[step] = reward;
	
            total_reward += reward;
            discounted_reward += discount * reward;
			
            discount *= gamma;            
			
            episode++;
            statistics.ep_stats.resize(episode + 1);
            statistics.ep_stats[episode].total_reward = 0.0;
            statistics.ep_stats[episode].discounted_reward = 0.0;
            statistics.ep_stats[episode].steps = 0;
            discount = 1.0;
			printf ("# episode %d complete Step = %d\n", episode,current_time);
            
			environment->Reset();

            action_ok = true;
			current_time = 0;
            if (n_episodes >= 0 && episode >= n_episodes) {
				// logmsg (" Breaking after %d episodes,  %d steps\n", episode, step);
                break;
            }
            step++;
        }
        Vector state = environment->getState();
		int action = algorithm->Act(state);		
        real reward = environment->getReward();
        
		statistics.reward.resize(step + 1);
        statistics.reward[step] = reward;
        statistics.ep_stats[episode].steps++;
        statistics.ep_stats[episode].total_reward += reward;
        statistics.ep_stats[episode].discounted_reward += discount * reward;
        total_reward += reward;
        discounted_reward += discount * reward;
		
        discount *= gamma;
		
        action_ok = environment->Act(action);
        current_time++;
		
    }
	//    printf(" %f %f # RUN_REWARD\n", total_reward, discounted_reward);
	
    if ((int) statistics.ep_stats.size() != n_episodes) {
        statistics.ep_stats.resize(statistics.ep_stats.size() - 1);
    }
	// printf ("# Exiting after %d episodes, %d steps (%d %d)\n",
	//            episode, n_steps,
	//            (int) statistics.ep_stats.size(),
	//            (int) statistics.reward.size());
	//		
    return statistics;
}


#endif
