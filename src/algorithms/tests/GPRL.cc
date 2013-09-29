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
#include "SparseGaussianProcess.h"
//#include "SparseGreedyGaussianProcess.h"
#include "GP_FittedLSTD.h"



// -- Usual includes -- //
#include <cstring>
#include <getopt.h>
/** Options */
struct Options
{
    real gamma;						///< discount factor
    const char* environment_name;	///< environment name
    RandomNumberGenerator& rng;		///< random number generator
    int n_samples;					///< number of sampled pairs for LSTD training
	int n_training_e_steps;			///< Max length of training trajectories
    int n_training_episodes;		///< number of training trajectories
	int n_testing_e_steps;			///< Max length of testing trajectories.
    int n_testing_episodes;			///< number of testing trajectories
    int grid;						///< grid size for features
    real scale;						///< scale of RBF basis
	real randomness;				///< environment randomness
    Options(RandomNumberGenerator& rng_) :
	gamma(0.999),
	environment_name("MountainCar"),
	rng(rng_), 
	n_samples(3000),
	n_training_e_steps(40),
	n_training_episodes(1000),
	n_testing_e_steps(1000),
	n_testing_episodes(1000),
	grid(4),
	scale(1.0),
	randomness(0.0)
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
        logmsg("n_samples: %d\n", n_samples);
        logmsg("n_training_e_steps: %d\n",  n_training_e_steps);
		logmsg("n_training_episodes: %d\n", n_training_episodes);
		logmsg("n_testing_e_steps: %d\n", n_testing_e_steps);
		logmsg("n_testing_episodes: %d\n", n_testing_episodes);
        logmsg("grid %d\n", grid);
        logmsg("scale %f\n", scale);
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

void OfflineGPRL(real gamma,
				 int n_steps, 
				 int n_episodes,
				 GP_FittedLSTD<Vector,int>* algorithm,
				 std::vector<std::vector<SparseGaussianProcess*> > EnvPrediction,
				 ContinuousStateEnvironment* environment);

Statistics OnlineGPRL(real gamma,
					  int n_steps, 
					  int n_episodes,
					  GP_FittedLSTD<Vector,int>* algorithm,
					  std::vector<std::vector<SparseGaussianProcess*> > EnvPrediction,
					  ContinuousStateEnvironment* environment);

void EvaluatePrediction(int N,
						std::vector<std::vector<SparseGaussianProcess*> > algorithm,
						ContinuousStateEnvironment* environment);

Statistics EvaluateAlgorithm (real gamma,
							  int episode_steps,
                              int n_episodes,
                              GP_FittedLSTD<Vector,int>* algorithm,
                              ContinuousStateEnvironment* environment);

static const char* const help_text = "Usage: online_algorithms [options] algorithm environment\n\
\nOptions:\n\
	--environment:			{MountainCar, Pendulum, Puddle, Bicycle, CartPole, Acrobot}\n\
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
\n\
* denotes default parameter\n\
\n";


int main(int argc, char* argv[])
{
	ulong seed = time(NULL);
	char* seed_filename = 0;
	const char* test_mode = "offline";
	
	bool online_test = false;
		
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
				{"n_samples", required_argument, 0, 0},				//2
				{"n_training_e_steps", required_argument, 0, 0},	//3
				{"n_training_episodes", required_argument, 0, 0},	//4
				{"n_testing_e_steps", required_argument, 0, 0},		//5
				{"n_testing_episodes", required_argument, 0, 0},	//6
				{"grid", required_argument, 0, 0},					//7
				{"scale", required_argument, 0, 0},					//8
				{"online", no_argument, 0, 0},						//9
				{"seed", required_argument, 0, 0},					//10
				{"seed_file", required_argument, 0, 0},				//11
				{"randomness", required_argument, 0, 0},			//12
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
						case 2: options.n_samples = atoi(optarg); break;
						case 3: options.n_training_e_steps = atoi(optarg); break;
						case 4: options.n_training_episodes = atoi(optarg); break;
						case 5: options.n_testing_e_steps = atoi(optarg); break;
						case 6: options.n_testing_episodes = atoi(optarg); break;
						case 7: options.grid = atoi(optarg); break;
						case 8: options.scale = atof(optarg); break;
						case 9: online_test = true; break;
						case 10: seed = atoi(optarg); break;
						case 11: seed_filename = optarg; break;
						case 12: options.randomness = atof(optarg); break;
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
	
	std::cout << "Basis functions creation" <<std::endl;
	RBFBasisSet* RBFs = NULL;
	std::cout << "Creating Radial basis functions..." << std::endl;
	EvenGrid Discretisation(S_L, S_U, options.grid);
	
	RBFs = new RBFBasisSet(Discretisation, options.scale);
		
	int m = RBFs->size() + 1; // redefinition of the input dimensions (size of the created basis functions plus a dummy state)
	std::cout << "# Number of basis functions: " << m << std::endl;

	
	std::cout << "# Creating " << n_actions*n_states << " Gaussian Prediction environement models..." << std::endl;
	
	//bayesian multivariate regression model for the system transition model (one for each action & dimension)
	std::vector<std::vector<SparseGaussianProcess*> > GaussianPrediction(n_actions, std::vector<SparseGaussianProcess*>(n_states));
	
	///Hyperparameters.
	real noise_variance = 1e-10;
	Vector scale_length = abs(S_L - S_U) / 10;
	Vector scale_length_dic = abs(S_L - S_U) / 20;
	real sig_val = 1.0;
	for ( int i = 0; i < n_actions; ++i) {
		for( int j = 0; j < n_states; ++j) {
			GaussianPrediction[i][j] = new SparseGaussianProcess( noise_variance, scale_length, sig_val, 0.1, scale_length_dic);
		}
	}
	std::cout << "Creation of the Gaussian Prediction environement models completed..." << std::endl;
	
	std::cout << "Fitted LSTD Initialization..." << std::endl;
	GP_FittedLSTD<Vector,int> *algorithm = new GP_FittedLSTD<Vector,int>(options.gamma, options.n_samples, environment, GaussianPrediction, RBFs);

	int episodes;
	Statistics run_statistics;
	if(online_test) {	
		test_mode ="online";
		episodes = options.n_training_episodes;
		run_statistics = OnlineGPRL(options.gamma,
									options.n_testing_e_steps, 
									options.n_testing_episodes,
									algorithm,
									GaussianPrediction,
									environment);
	} else {
		episodes = options.n_testing_episodes;
		OfflineGPRL(options.gamma,
					options.n_training_e_steps,
					options.n_training_episodes,
					algorithm,
					GaussianPrediction,
					environment);
		
		run_statistics = EvaluateAlgorithm(options.gamma,
										   options.n_testing_e_steps,
										   options.n_testing_episodes,
										   algorithm, 
										   environment);
	}
	
	Vector Stats(episodes);
	Vector Statr(episodes);
	real train_steps = 0.0;
	for( uint i = 0; i < run_statistics.ep_stats.size(); ++i) {
		Stats(i) = run_statistics.ep_stats[i].steps;
		Statr(i) = run_statistics.ep_stats[i].total_reward;
		train_steps += run_statistics.ep_stats[i].steps;
	}
	printf("Mean number of steps = %f\n", train_steps / episodes);
			
	char buffer[100];
	sprintf (buffer, "%s_GPRL_RESULTS_STEPS_%s",test_mode, options.environment_name);
	FILE *output	= fopen(buffer,"w");
	if(output!=NULL) { 
		Stats.print(output);
	}
	fclose(output);
	sprintf (buffer, "%s_GPRL_RESULTS_REWARDS_%s",test_mode, options.environment_name);
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

void OfflineGPRL(real gamma,
				 int n_steps, 
				 int n_episodes,
				 GP_FittedLSTD<Vector,int>* algorithm,
				 std::vector<std::vector<SparseGaussianProcess*> > EnvPrediction,
				 ContinuousStateEnvironment* environment)
{
	std::cout << "#Offline training ..." << environment->Name() << std::endl;

	Vector state, next_state;
	int action;
	real reward;
	int n_actions = environment->getNActions();
	int n_states  = environment->getNStates();
	std::vector< std::vector<Vector> > Input(n_actions);
	std::vector< std::vector<real> > Output(n_actions*n_states);

	for(int episode = 0; episode < n_episodes; ++episode) {
		environment->Reset();
		state = environment->getState();

		action = (int) floor(urandom(0.0, (real) environment->getNActions()));
		int step		= 0;
		bool action_ok	= true;
		
		while(action_ok && step < n_steps) {
			action_ok	= environment->Act(action);
			
			reward		= environment->getReward();
			next_state	= environment->getState();
			
			Input[action].push_back(state);
			for(int dim = 0; dim<n_states; ++dim) {
				Output[(action*n_states + dim)].push_back(next_state(dim));
			}
			action = (int)floor(urandom(0.0, (real) environment->getNActions()));
			state = next_state;
			step++;
		}
	}
	/// In this point we update the gaussian processes
	for(int a = 0; a<n_actions; ++a) {
		for(int dim = 0; dim<n_states; ++dim) {
			EnvPrediction[a][dim]->Observe(Input[a], Output[a*n_states + dim]);
		}
	}
	algorithm->Update();
//	EvaluatePrediction(5000,
//						  EnvPrediction,
//						  environment);
	printf("The offline training is completed\n");
}

Statistics OnlineGPRL(real gamma, 
					  int n_steps, 
					  int n_episodes,
					  GP_FittedLSTD<Vector,int>* algorithm,
					  std::vector<std::vector<SparseGaussianProcess*> > EnvPrediction,
					  ContinuousStateEnvironment* environment)
{
	std::cout << "#Online training ..." << environment->Name() << std::endl;
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
			
			Input[action].push_back(state);
			for(int dim = 0; dim<n_states; ++dim) {
				Output[(action*n_states + dim)].push_back(next_state(dim));
			}
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
		/// In this point we update the gaussian processes
		for(int a = 0; a<n_actions; ++a) {
			if(Input[a].size() > 0) {
				for(int dim = 0; dim<n_states; ++dim) {
					EnvPrediction[a][dim]->AddObservation(Input[a], Output[a*n_states + dim]);
					Output[a*n_states + dim].clear();
				}
				Input[a].clear();
			}
		}
		algorithm->Update();
	}
	printf("The online training is completed\n");

	return statistics;
}

void EvaluatePrediction(int N,
						   std::vector<std::vector<SparseGaussianProcess*> > algorithm,
						   ContinuousStateEnvironment* environment)
{
	Vector S_L	= environment->StateLowerBound();
	Vector S_U	= environment->StateUpperBound();
	int n_actions	= environment->getNActions();
	int n_states	= environment->getNStates();
	
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
	
	for(int a=0; a<n_actions; ++a) {
		for(int dim=0; dim<n_states; ++dim) {
			n = sprintf (buffer, "Output_samples_action_%d_Dimension_%d", a,dim);
			np = sprintf(bufferp, "Predicted_Output_Action_%d_Dimension_%d",a,dim);
			FILE *output	= fopen(buffer,"w");
			FILE *outputp   = fopen(bufferp,"w");
			if(output!=NULL) {
				for( uint i = 0; i < states.size(); ++i) {
					environment->Reset();
					environment->setState(states[i]);
					environment->Act(a);
					environment->getState().print(output);
					fprintf(outputp,"%f \n", algorithm[a][dim]->GeneratePrediction(states[i]));
				}
			}
			fclose(output);
			fclose(outputp);
		}
	}
}

/*** Evaluate an algorithm
episode_steps: maximum number of steps per episode. If negative, then ignore
n_steps: maximun number of total steps. If negative, then ignore.
n_episodes: maximum number of episodes. Cannot be negative.
*/
Statistics EvaluateAlgorithm (real gamma,
							  int episode_steps,
                              int n_episodes,
                              GP_FittedLSTD<Vector,int>* algorithm,
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