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
#include "BayesianMultivariate.h"
#include "GaussianProcess.h"
#include "GP_FittedLSTD.h"



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

void EnvironmentPrediction(int N,
						   std::vector<Vector> Samples,
						   int action,
						   int dim,
						   GaussianProcess* algorithm,
						   ContinuousStateEnvironment* environment);

Statistics EvaluateAlgorithm (int episode_steps,
                              int n_episodes,
                              GP_FittedLSTD<Vector,int>* algorithm,
                              ContinuousStateEnvironment* environment,
                              real gamma);

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
--environment:		{PuddleWorld, MountainCar, Pendulum, Linear}\n\
--n_states:			number of states (usually there is no need to specify it)\n\
--n_actions:		number of actions (usually there is no need to specify it)\n\
--gamma:			reward discounting in [0,1] (* 0.99)\n\
--randomness:		environment randomness (* 0.0)\n\
--n_runs:			maximum number of runs (* 1)\n\
--n_train_episodes:	maximum number of train episodes (ignored if < 0)\n\
--n_test_episodes:  maximum number of test episodes (ignored if < 0)\n\
--episode_steps:    maximum number of steps in each episode (ignored if <0)\n\
--n_steps:			maximum number of total steps\n\
--rbf:				use of basis functions or not (* 0 or 1)\n\
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
	int n_states			= 2;
	real gamma				= 0.999;
	real randomness			= 0.0;
	real epsilon			= 0;
	uint n_runs				= 1;
    uint grids				= 4;
	uint n_test_steps       = 1000;
	uint n_test_episodes    = 1000;
	int n_samples			= 100;
	int sampling			= 1;
	real rbf_scale			= 1.0;
	//Value iteration threshold
	real threshold			= 0.001;
	
	
	const char * environment_name = "Pendulum";
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
				{"n_runs", required_argument, 0, 0},			//3
				{"n_test_steps", required_argument, 0, 0},		//4
				{"n_test_episodes", required_argument, 0, 0},   //5
				{"n_samples", required_argument, 0, 0},			//6
                {"epsilon", required_argument, 0, 0},			//7
                {"environment_name", required_argument, 0, 0},	//8
                {"grids", required_argument, 0, 0},				//9
                {"randomness", required_argument, 0, 0},		//10
				{"sampling", required_argument, 0, 0},			//11
				{"rbf_scale", required_argument, 0, 0},			//12
				{"threshold",required_argument, 0, 0},			//13
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
						case 0: n_states = atoi(optarg); break;
						case 1: n_actions = atoi(optarg); break;
						case 2: gamma = atof(optarg); break;
						case 3: n_runs = atoi(optarg); break;
						case 4: n_test_steps = atoi(optarg); break;
						case 5: n_test_episodes = atoi(optarg); break;
						case 6: n_samples = atoi(optarg); break;
						case 7: epsilon = atof(optarg); break; 
						case 8: environment_name = optarg; break;
						case 9: grids = atoi(optarg); break;
						case 10: randomness = atof(optarg); break;
						case 11: sampling = atoi(optarg); break;
						case 12: rbf_scale = atoi(optarg); break;
						case 13: threshold = atof(optarg); break;
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
	assert(n_states > 0);
    assert(n_actions > 0);
    assert(gamma >= 0 && gamma <= 1);
    assert(randomness >= 0 && randomness <= 1);
    assert(n_runs > 0);
	assert(n_test_steps > 0);
	assert(n_test_episodes > 0);
    assert(n_samples > 0);
	assert(rbf_scale > 0);
    assert(grids > 0);
		
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
	
    std::cout << "Starting GPRL" << std::endl;
    
   	
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
	}
	else if (!strcmp(environment_name, "Acrobot")) {
		environment = new Acrobot();
	}
	else if (!strcmp(environment_name, "CartPole")) {
		environment = new CartPole();
	}
	else if (!strcmp(environment_name, "Bike")) {
		environment = new Bike();
	}
	else if(!strcmp(environment_name,"Linear")) {
		environment = new LinearDynamicQuadratic();
	}
	else {
		fprintf(stderr, "Unknown environment %s \n", environment_name);
	}
	
	std::cout << environment_name << " environment creation completed" << std::endl;
	// making sure the number of states & actions is correct
	n_states	= environment->getNStates();
	n_actions	= environment->getNActions();
	Vector S_L	= environment->StateLowerBound();
	Vector S_U	= environment->StateUpperBound();
	
	std::cout <<  "Creating environment: " << environment_name
	<< " with " << n_states << " states and , "
	<< n_actions << " actions.\n";
	
	std::cout << "Basis functions creation" <<std::endl;
	RBFBasisSet* RBFs = NULL;
	std::cout << "Creating Radial basis functions..." << std::endl;
	EvenGrid Discretisation(S_L, S_U, grids);
	
	RBFs = new RBFBasisSet(Discretisation,rbf_scale);
		
	int m = RBFs->size() + 1; // redefinition of the input dimensions (size of the created basis functions plus a dummy state)
	std::cout << "# Number of basis functions: " << m << std::endl;

	
	std::cout << "# Creating " << n_actions*n_states << " Gaussian Prediction environement models..." << std::endl;
	
	//bayesian multivariate regression model for the system transition model
	//one for each action
	std::vector<std::vector<GaussianProcess*> > GaussianPrediction(n_actions, std::vector<GaussianProcess*>(n_states));
//	GaussianPrediction.resize(n_actions);
	
	///Hyperparameters.
	real noise_variance = 0.1;
	Vector scale_length = abs(S_L - S_U)/20;
	real sig_val = 1.0;
	for ( int i = 0; i < n_actions; ++i) {
		for( int j = 0; j < n_states; ++j) {
			GaussianPrediction[i][j] = new GaussianProcess(noise_variance, scale_length, sig_val);
		}
	}
	std::cout << "Creation of the Gaussian Prediction environement models completed..." << std::endl;
	
	std::cout << "Fitted LSTD Initialization..." << std::endl;
	GP_FittedLSTD<Vector,int> *algorithm = new GP_FittedLSTD<Vector,int>(gamma, 5000, environment, GaussianPrediction, RBFs);

	Matrix Stats(n_runs, n_test_episodes);
	Matrix Statr(n_runs, n_test_episodes);
	
	for ( uint run = 0; run < n_runs; ++run) {
		std::cout << "# Run: " << run << std::endl;
		
		std::cout << "#Collecting " << n_samples << " uniformly random samples for system dynamics identification" << std::endl;
		std::vector<Vector> Samples;
		for(int i=0; i<n_samples; ++i) {
			Vector state = urandom(S_L, S_U);
			Samples.push_back(state);
		}
		
		for(int i = 0; i < n_actions; ++i) {
			for(int j = 0; j < n_states; ++j) {
				EnvironmentPrediction(n_samples, Samples, i, j, GaussianPrediction[i][j], environment);
			}
		}
		
		algorithm->Update();

		epsilon = 0.0;
		Statistics run_statistics = EvaluateAlgorithm(n_test_steps, 
													  n_test_episodes,
													  algorithm, 
													  environment,
													  gamma);
		real train_steps = 0.0;
		for( uint i = 0; i < run_statistics.ep_stats.size(); ++i) {
			Stats(run,i) = run_statistics.ep_stats[i].steps;
			Statr(run,i) = run_statistics.ep_stats[i].total_reward;
			train_steps += run_statistics.ep_stats[i].steps;
		}
		printf("Mean number of steps = %f\n", train_steps / n_test_episodes);
		//for (uint i=0; i<run_statistics.ep_stats.size(); ++i) {
		//		
		//			statistics.ep_stats[i].total_reward += run_statistics.ep_stats[i].total_reward;
		//			statistics.ep_stats[i].discounted_reward += run_statistics.ep_stats[i].discounted_reward;
		//			statistics.ep_stats[i].steps += run_statistics.ep_stats[i].steps;
		//			statistics.ep_stats[i].mse += run_statistics.ep_stats[i].mse;
		//			statistics.ep_stats[i].n_runs++;
		//		}
		//
		//		for (uint i=0; i<run_statistics.reward.size(); ++i) {
		//			statistics.reward[i] += run_statistics.reward[i];
		//			statistics.n_runs[i]++;
		//		}
		
		algorithm->Reset();
	}
	
	char buffer[100];
	sprintf (buffer, "GPRL_RESULTS_STEPS_%s->(#Samples = %d)",environment_name,n_samples);
	FILE *output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Stats.print(output);
	}
	fclose(output);
	sprintf (buffer, "GPRL_RESULTS_REWARDS_%s->(#Samples = %d)",environment_name,n_samples);
	output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Statr.print(output);
	}
	fclose(output);
	
	for( int i = 0; i<n_actions; ++i){
		for(int j = 0; j<n_states; ++j) {
			delete GaussianPrediction[i][j];
		}
	}
	delete algorithm;
	delete environment;
	
	std::cout << "Done" << std::endl;
	
	return 0;
}

void EnvironmentPrediction(int N,
						   std::vector<Vector> Samples, 
						   int action,
						   int dim,
						   GaussianProcess* algorithm,
						   ContinuousStateEnvironment* environment)
{
	std::cout << "#Environment System dynamics identification for Action => " << action << " and Dimension => " << dim << std::endl;
	Vector state;
	Vector next_state;
	Vector S_L	= environment->StateLowerBound();
	Vector S_U	= environment->StateUpperBound();
	
	Matrix X(N, environment->getNStates()); 
	Vector Y(N); 
	for(int i=0; i<N; ++i) {
		X.setRow(i, Samples[i]);
		environment->Reset();
		environment->setState(Samples[i]);
		// We follow the action
		environment->Act(action);

		next_state	= environment->getState();
		Y(i) = next_state(dim);
	}
	algorithm->Observe(X,Y);
	
	///Prediction model testing.
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
	char bufferp[100];
	int n, np;

	n = sprintf (buffer, "Output_samples_action_%d_Dimension_%d", action,dim);
	np = sprintf(bufferp, "Predicted_Output_Action_%d_Dimension_%d",action,dim);
	FILE *output	= fopen(buffer,"w");
	FILE *outputp   = fopen(bufferp,"w");
	if(output!=NULL) {
		for( uint i = 0; i < states.size(); ++i) {
			environment->Reset();
			environment->setState(states[i]);
			environment->Act(action);
			environment->getState().print(output);
			fprintf(outputp,"%f \n", algorithm->GeneratePrediction(states[i]));
		}
	}
	fclose(output);
	fclose(outputp);

	printf("Model training end\n");
}

/*** Evaluate an algorithm

episode_steps: maximum number of steps per episode. If negative, then ignore
n_steps: maximun number of total steps. If negative, then ignore.
n_episodes: maximum number of episodes. Cannot be negative.
*/

Statistics EvaluateAlgorithm (int episode_steps,
                              int n_episodes,
                              GP_FittedLSTD<Vector,int>* algorithm,
                              ContinuousStateEnvironment* environment,
                              real gamma)
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
			//  if (episode >= 0) {
			//                statistics.ep_stats[episode].steps++;
			//                statistics.ep_stats[episode].total_reward += reward;
			//				statistics.ep_stats[episode].discounted_reward += discount * reward;
			//            }
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
			environment->getState(); //.print(stdout);
			
            action_ok = true;
			current_time = 0;
            if (n_episodes >= 0 && episode >= n_episodes) {
				// logmsg (" Breaking after %d episodes,  %d steps\n", episode, step);
                break;
            }
            step++;
        }
		
        Vector state = environment->getState();
		//        state.print(stdout);
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
		
		//		printf("action = %d\n",action);
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