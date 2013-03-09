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

#include "OnlineLSPI.h"
#include "MersenneTwister.h"
#include "Grid.h"
#include "BasisSet.h"
#include "Rollout.h"
#include "Pendulum.h"
#include "PuddleWorld.h"
#include "CartPole.h"
#include "Pendulum.h"
#include "Acrobot.h"
#include "MountainCar.h"
#include "MersenneTwister.h"
#include "RandomPolicy.h"
#include "Random.h"
#include "ContinuousPolicy.h"
#include "Vector.h"
#include <getopt.h>
#include <cstring>

// -- Randomness -- //
#include "RandomNumberFile.h"

struct PerformanceStatistics
{
    real total_reward;
    real discounted_reward;
    real run_time;
    PerformanceStatistics()
	: total_reward(0),
	discounted_reward(0),
	run_time(0)
    {}
};

struct AveragePerformanceStatistics : PerformanceStatistics
{
    int N;
    AveragePerformanceStatistics()
	: N(0)
    {
    }
    void Observe(PerformanceStatistics& x)
    {
        N++;
        real rN = (real) N;
        real irN = 1 / rN;
        real gamma = (rN - 1) * irN;
		
        total_reward = gamma * total_reward + (1 - gamma) * x.total_reward;
        discounted_reward = gamma * discounted_reward + (1 - gamma) * x.discounted_reward;
        run_time = gamma * run_time + (1 - gamma) * x.run_time;
		
    }
    void Show()
    {
        printf ("%f %f %f # reward, discounted, time\n", 
                total_reward,
                discounted_reward,
                run_time);
    }
};

PerformanceStatistics Evaluate(Environment<Vector, int>* environment,
                               AbstractPolicy<Vector, int>& policy,
                               real gamma, int T);

static const char* const help_text = "Usage: rsapi [options]\n\
\nOptions:\n\
--environment: {PuddleWorld, Pendulum, Mountain}\n\
--gamma:       reward discounting in [0,1]\n\
--max_iteration:      maximum number of policy iterations (lspi)\n\
--horizon:     rollout horizon\n\
--delta:       stopping threshold\n\
--grids:       RBF factor\n\
--n_rollouts:  number of rollouts\n\
--n_runs:		   number of runs\n\
--eval_iter:   number of evaluations\n\
--algorithm: {1:LSTDQ, 2:LSTD1_OPT}\n\
\n";

int main(int argc, char* argv[])
{
	//Create a new environment.
	Environment<Vector,int>* environment;
	
	int seed = time(NULL);
	
	int max_iteration	= 20;
	real gamma	       	= 0.999;
	int episodes            = 1000;
	uint episode_steps	= 1000;
	int grids		= 5;
	int algorithm		= 1;
	real delta		= 0.0001;
	char* environment_name  = NULL;
	int eval_iter		= 1000;
	int n_runs		= 100;
	real scale		= 1.0; 
	
	{
		//options
		int c;
		int digit_optind = 0;
		while(1){
			int this_option_optind = optind ? optind : 1;
			int option_index = 0;
			static struct option long_options[] = { 
				{"max_iteration",required_argument, 0, 0}, //0
				{"gamma",required_argument, 0, 0}, //1
				{"episodes",required_argument, 0 ,0}, //2
				{"episode_steps", required_argument, 0 , 0}, // 3
				{"environment", required_argument, 0, 0}, //4
				{"delta", required_argument, 0, 0}, //5
				{"grids", required_argument, 0, 0}, //6
				{"n_runs",required_argument, 0, 0}, //7
				{"eval_iter", required_argument, 0, 0}, //8
				{"algorithm", required_argument, 0, 0}, //9
				{"scale", required_argument, 0, 0}, //10
				{"seed", required_argument, 0, 0}, //11
				{0, 0, 0, 0}
			};
			c = getopt_long(argc, argv, "", long_options, &option_index);
			if(c == -1)
			  break;
			
			switch (c){
			case 0:
#if 0
	       		    printf ("option %s (%d)", long_options[option_index].name, option_index);
	       		    if (optarg)
	       		       	printf (" with arg %s", optarg);
	       		    printf ("\n");
#endif	
	       		    switch (option_index) {
			    case 0: max_iteration = atoi(optarg); break;
			    case 1: gamma = atof(optarg); break;
			    case 2: episodes = atoi(optarg); break;
			    case 3: episode_steps = atoi(optarg); break;
			    case 4: environment_name = optarg; break;
			    case 5: delta = atof(optarg); break;
			    case 6: grids = atoi(optarg); break;
			    case 7: n_runs = atoi(optarg); break;
			    case 8: eval_iter = atoi(optarg); break;
			    case 9: algorithm = atoi(optarg); break;
			    case 10: scale = atof(optarg); break;
			    case 11: seed = atoi(optarg); break;
			    default:
			        fprintf (stderr, "Invalid options\n");
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
	
	RandomNumberGenerator* rng;
	RandomNumberGenerator* environment_rng;
	
	MersenneTwisterRNG mersenne_twister_env;
	environment_rng = (RandomNumberGenerator*) &mersenne_twister_env;
	
	MersenneTwisterRNG mersenne_twister;
	rng = (RandomNumberGenerator*) &mersenne_twister;
	
	
	if (!environment_name) {
	  fprintf(stderr, "Must specify environment\n");
	  exit(-1);
	}
	if (!strcmp(environment_name, "PuddleWorld")) {
	  environment = new PuddleWorld();
	} else if (!strcmp(environment_name, "Pendulum")) {
	  environment = new Pendulum(); //2.0, 8.0, 0.5, 9.8, 0.0);    
	} else if (!strcmp(environment_name, "CartPole")) {
	  environment = new CartPole(); //2.0, 8.0, 0.5, 9.8, 0.0);    
	} else if (!strcmp(environment_name, "Acrobot")) {
	  environment = new Acrobot(); //2.0, 8.0, 0.5, 9.8, 0.0);    
	} else if (!strcmp(environment_name, "Mountain")) {
	  environment = new MountainCar();
	} else {
	  fprintf(stderr, "Invalid environment name %s\n", environment_name);
	  exit(-1);
	}
		
	int state_dimension = environment->getNStates();
	int n_actions = environment->getNActions();
	Vector S_L = environment->StateLowerBound();
	Vector S_U = environment->StateUpperBound();
	
	printf("# State dimension: %d\n", state_dimension);
	printf("# S_L: "); S_L.print(stdout);
	printf("# S_U: "); S_U.print(stdout);
//      S_U[0] = (3.0*M_PI)/4.0; //4;
//	S_U[1] = 1.5;//10;
//	S_L[0] = (-3.0*M_PI)/4.0;//-4;
//	S_L[1] = -1.5;//10;
//      Vector D = Vector::Unity(2);
//	EvenGrid Discretisation(S_L, S_U, D, grids);
	EvenGrid Discretisation(S_L, S_U, grids);

	srand48(34987235);
	srand(34987235);
	setRandomSeed(34987235);
	environment_rng->manualSeed(228240153);
	rng->manualSeed(1361690241);	
	
	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;
	
	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), rng);
	
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout
	= new Rollout<Vector, int, AbstractPolicy<Vector, int> >( environment->getState(), policy, environment, gamma, true);
	
	rollout->Sampling(episodes, episode_steps);
//	printf("Total number of collected samples -> %d\n",rollout->getNSamples());

	
	RBFBasisSet* RBFs = new RBFBasisSet(Discretisation,scale);
	OnlineLSPI* lspi = new OnlineLSPI(gamma, delta, state_dimension, n_actions, max_iteration, RBFs);
	
	Vector state;
	Vector next_state;
	int action;
	int next_action;
	bool action_ok;
	real reward;
	real steps;
	real total_reward;

	FixedContinuousPolicy policy_lspi = lspi->ReturnPolicy();
	
	real epsilon; //e-greedy
	
	Matrix Stats(n_runs, episodes);
	Matrix Statr(n_runs, episodes);
	
	for(int run = 0; run < n_runs; ++run) {
	  std::cout << "# Run: " << run << std::endl;
	  epsilon = 1.0;
	  for (int ep=0; ep<episodes; ++ep) {
	    environment->Reset();
	    state = environment->getState();
			
	    action = policy_lspi.SelectAction(state);
			
	    steps = 0;
	    total_reward = 0;
	    for(uint t = 0; t<episode_steps; ++t) {
	      action_ok	= environment->Act(action);
	      reward		= environment->getReward();
			
	      next_state = environment->getState();
	      if(urandom() < epsilon)
		next_action = (int)urandom(0, environment->getNActions());
	      else
		next_action = policy_lspi.SelectAction(next_state);
				
	      lspi->LSTDQ(state, action, reward, next_state, next_action, environment->getEndsim(),false);
	      steps++;
	      total_reward += reward;
	      epsilon *= 0.9968;
				
	      state		= next_state;
	      action		= next_action;
	      if(!action_ok) {
		break;
	      }
	    }
	    policy_lspi = lspi->ReturnPolicy();

	    lspi->Update();
	    printf("Episode %d completed (Steps = %f)\n",ep, steps);
	    Stats(run,ep) = steps;
	    Statr(run,ep) = total_reward;
	  }
	  char buffer[100];
	  sprintf (buffer, "LSPI_ONLINE_RESULTS_STEPS_%s",environment_name);
	  FILE *output	= fopen(buffer,"w");
	  if(output!=NULL) {
	    Stats.print(output);
	  }
	  fclose(output);
	  sprintf (buffer, "LSPI_ONLINE_RESULTS_REWARDS_%s",environment_name);
	  output = fopen(buffer,"w");
	  if(output!=NULL) {
	    Statr.print(output);
	  }
	  fclose(output);
	  lspi->Reset();
	}
	
	delete policy;
			
	delete lspi;
	delete rollout;
			
	delete RBFs;
	
	delete environment;
	
}

PerformanceStatistics Evaluate(Environment<Vector, int>* environment,
							   AbstractPolicy<Vector, int>& policy,
							   real gamma, int T)
{
	PerformanceStatistics statistics;
	statistics.total_reward = 0;
	statistics.discounted_reward = 0;
	statistics.run_time = 0;
	int t;
	real discount = 1;
	environment->Reset();
	policy.Reset();
//	printf("New episode\n");
	for (t=0; t < T; ++t, ++statistics.run_time) {
	  Vector state = environment->getState();
//		state.print(stdout);
	  real reward = environment->getReward();
		//        printf("%d %f # r_t\n", t, reward);
	  if (t > 0) {
            statistics.total_reward += reward;
	    statistics.discounted_reward += reward*discount;
            discount *=gamma;
	  }
	  policy.setState(state);
	  int action = policy.SelectAction();
	  bool action_ok = environment->Act(action);
	  if(!action_ok) {
            real reward = environment->getReward();
	    statistics.total_reward += reward;
            statistics.discounted_reward += reward*discount;
	    break;
	  }
	}
//	printf("Run time = %f\n",statistics.run_time);
	return statistics;
}

#endif


