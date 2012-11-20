/*
 *  lspi.cc
 *  
 *
 *  Created by Νικόλαος Τζιωρτζιωτης on 18/11/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef MAKE_MAIN
#include "LSPI.h"
#include "Grid.h"
#include "BasisSet.h"
#include "Rollout.h"
#include "Pendulum.h"
#include "PuddleWorld.h"
#include "Pendulum.h"
#include "RandomPolicy.h"
#include "Random.h"
#include "MersenneTwister.h"
#include "ContinuousPolicy.h"
#include "Vector.h"
#include <getopt.h>
#include <cstring>

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
--environment: {PuddleWorld, Pendulum}\n\
--gamma:       reward discounting in [0,1]\n\
--max_iteration:      maximum number of policy iterations\n\
--horizon:     rollout horizon\n\
--delta:       stopping threshold\n\
--grids:       RBF factor\n\
--n_rollouts:  number of rollouts\n\
\n";

int main(int argc, char* argv[])
{
	//Create a new environment.
	Environment<Vector,int>* environment;
	
	MersenneTwisterRNG rng;

	int max_iteration = 100;
	real gamma = 0.999;
	int n_rollouts = 100;
	int horizon = 50;
	int grids = 10;
	real delta = 0.01;
	char* environment_name = NULL;

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
				{"n_rollouts",required_argument, 0 ,0}, //2
				{"horizon", required_argument, 0 , 0}, // 3
				{"environment", required_argument, 0, 0}, //4
				{"delta", required_argument, 0, 0}, //5
				{"grids", required_argument, 0, 0}, //6
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
						case 2: n_rollouts = atoi(optarg); break;
						case 3: horizon = atoi(optarg); break;
						case 4: environment_name = optarg; break;
						case 5: delta = atof(optarg); break;
						case 6: grids = atoi(optarg); break;
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
	
	if (!environment_name) {
        fprintf(stderr, "Must specify environment\n");
        exit(-1);
    }
    if (!strcmp(environment_name, "PuddleWorld")) {
        environment = new PuddleWorld();
    } else if (!strcmp(environment_name, "Pendulum")) {
        environment = new Pendulum();    
    } else {
        fprintf(stderr, "Invalid environment name %s\n", environment_name);
        exit(-1);
    }
	
	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;

	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), &rng);

	int state_dimension = environment->getNStates();
	int n_actions = environment->getNActions();
	Vector S_L = environment->StateLowerBound();
    Vector S_U = environment->StateUpperBound();
    
    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
	
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout
        = new Rollout<Vector, int, AbstractPolicy<Vector, int> >(urandom(S_L,S_U), policy, environment, gamma, true);
	
    rollout->Sampling(n_rollouts, horizon);

	printf("Total number of collected samples -> %d\n",rollout->getNSamples());
	EvenGrid Discretisation(S_L, S_U, grids);
	
	for(int i = 0; i < rollout->getNRollouts(); ++i)
	{
		printf("Total number of collected samples -> %d\n",rollout->getNSamples(i));
	}
	
	RBFBasisSet* RBFs = new RBFBasisSet(Discretisation);
	LSPI* lspi = new LSPI(gamma, delta, state_dimension, n_actions, max_iteration, RBFs,rollout);
	lspi->PolicyIteration();
	
	AveragePerformanceStatistics statistics;
	for (int i=0; i<100; ++i) {
		PerformanceStatistics run_statistics = Evaluate(environment,
														lspi->ReturnPolicy(),
														gamma,
														1000);
		statistics.Observe(run_statistics);
	}
	statistics.Show();
	delete policy;
	delete RBFs;
	delete environment;
	delete lspi;
    delete rollout;
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
    for (t=0; t < T; ++t, ++statistics.run_time) {
		Vector state = environment->getState();
		real reward = environment->getReward();
		statistics.total_reward += reward;
		statistics.discounted_reward += reward*discount;
		discount *=gamma;
		policy.setState(state);
		int action = policy.SelectAction();
		bool action_ok = environment->Act(action);
		if(!action_ok) {
			break;
		}
	}
	return statistics;
}

#endif


