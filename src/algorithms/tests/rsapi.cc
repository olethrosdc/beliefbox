/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "RSAPI.h"
#include "Rollout.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "RandomPolicy.h"
#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "KNNClassifier.h"
#include "ClassifierPolicy.h"
#include <cstring>
#include <getopt.h>


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
                               AbstractPolicy<Vector, int>* policy,
                               real gamma, int T);


static const char* const help_text = "Usage: rsapi [options]\n\
\nOptions:\n\
    --environment: {MountainCar, Pendulum}\n\
    --n_states:    number of states\n\
    --gamma:       reward discounting in [0,1]\n\
    --n_iter:      maximum number of policy iterations\n\
    --horizon:     rollout horizon\n\
    --n_rollouts:  number of rollouts\n\
    --knn:         number of nearest neighbours in KNN classifier\n\
    --group:       use grouped action training\n\
    --resample:    resample states from discounted state distribution\n\
    --lipschitz:   use a Lipschitz assumption with constant L\n\
    --uniform:     uniformly sample all states in representative set\n\
    --upper_bound: use upper bound to sample states in set\n\
    --error_bound: sample uniformly until error bound is attained\n\
\n";

int main(int argc, char* argv[])
{
    
	// Create a new environment
	Environment<Vector, int>* environment;

	MersenneTwisterRNG rng;


	enum SamplingMethod 
	{
		UNIFORM,
		UPPER_BOUND,
		ERROR_BOUND
	};
    int n_states = 100;
    real gamma = 0.99;
    int n_iter=10;
    int n_rollouts = 100;
    int horizon = -1;
    int n_neighbours = 1;
    char* environment_name = NULL;
    bool group_training = false;
    bool resample = false;
    real delta = 0.1;
	real Lipschitz = -1;
	enum SamplingMethod sampling_method = UPPER_BOUND;

    {
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"n_states", required_argument, 0, 0}, //0
                {"gamma", required_argument, 0, 0}, //1
                {"n_iter", required_argument, 0, 0}, //2
                {"n_rollouts", required_argument, 0, 0}, //3
                {"horizon", required_argument, 0, 0}, //4
                {"environment", required_argument, 0, 0}, //5
                {"knn", required_argument, 0, 0}, //6
                {"group", no_argument, 0, 0}, //7
                {"resample", no_argument, 0, 0}, //8
                {"delta", required_argument, 0, 0}, //9
                {"lipschitz", required_argument, 0, 0}, //10
                {"uniform", no_argument, 0, 0}, //11
                {"upper_bound", no_argument, 0, 0}, //12
                {"error_bound", no_argument, 0, 0}, //13
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
                case 1: gamma = atof(optarg); break;
                case 2: n_iter = atoi(optarg); break;
                case 3: n_rollouts = atoi(optarg); break;
                case 4: horizon = atoi(optarg); break;
                case 5: environment_name = optarg; break;
                case 6: n_neighbours = atoi(optarg); break;
                case 7: group_training = true; break;
                case 8: resample = true; break;
                case 9: delta = atof(optarg); break;
                case 10: Lipschitz = atof(optarg); break;
				case 11: sampling_method = UNIFORM; break;
				case 12: sampling_method = UPPER_BOUND; break;
				case 13: sampling_method = ERROR_BOUND; break;
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
    if (!strcmp(environment_name, "MountainCar")) {
        environment = new MountainCar();
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
    Vector S_L = environment->StateLowerBound();
    Vector S_U = environment->StateUpperBound();
    
    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
    KNNClassifier* classifier = NULL;


    std::vector<Vector> state_vector(n_states);
    for (uint k=0; k<state_vector.size(); ++k) {
        Vector& state = state_vector[k];
        state.Resize(S_L.Size());
        for (int i=0; i<S_L.Size(); ++i) {
            state(i) = rng.uniform(S_L(i), S_U(i));
        }
    }

    //RSAPI rsapi(environment, &rng, gamma);
    for (int iter=0; iter<n_iter; ++iter) {
        AveragePerformanceStatistics statistics;
        for (int i=0; i<100; ++i) {
            PerformanceStatistics run_statistics = Evaluate(environment,
                                                            policy,
                                                            gamma,
                                                            1000);
            statistics.Observe(run_statistics);
        }
        statistics.Show();

        RSAPI rsapi(environment, &rng, gamma);
        rsapi.setPolicy(policy);
        for (uint k=0; k<state_vector.size(); ++k) {
            rsapi.AddState(state_vector[k]);
        }

		switch (sampling_method) {
		case UNIFORM:
			rsapi.SampleUniformly(n_rollouts, horizon);
			break;
		case ERROR_BOUND:
			rsapi.SampleToErrorBound(n_rollouts, horizon, delta);
			break;
		case UPPER_BOUND:
			rsapi.SampleUpperBound(n_rollouts, horizon, delta);
			break;
		}
		if (Lipschitz > 0) {
			rsapi.Bootstrap();
		}
        KNNClassifier* new_classifier = new KNNClassifier(state_dimension, environment->getNActions(), n_neighbours);
        int n_improved_actions = 0;
        if (group_training) {
            n_improved_actions = rsapi.GroupTrainClassifier(new_classifier, delta);
        } else {
            n_improved_actions = rsapi.TrainClassifier(new_classifier, delta);
        }
        printf ("# n: %d # improved actions\n", n_improved_actions);
        if (resample) {
            for (uint i=0; i<state_vector.size(); ++i) {
                if (rng.uniform() >= gamma) {
                    state_vector[i] =  rsapi.SampleStateFromPolicy();
                }
            }
        }
        delete policy;
        policy = new ClassifierPolicy(new_classifier);
        if (classifier) {
            delete classifier;
        }        
        classifier = new_classifier;
        fflush(stdout);
    }
    delete classifier;
    delete policy;
	delete environment;
}


PerformanceStatistics Evaluate(Environment<Vector, int>* environment,
                               AbstractPolicy<Vector, int>* policy,
                               real gamma, int T)
{
    PerformanceStatistics statistics;
    statistics.total_reward = 0;
    statistics.discounted_reward = 0;
    statistics.run_time = 0;

    real discount = 1;
    environment->Reset();
    for (int t=0; t < T; ++t, ++statistics.run_time) {
        Vector state = environment->getState();
        real reward = environment->getReward();
        statistics.total_reward += reward;
        statistics.discounted_reward += reward * discount;
        discount *= gamma;
        policy->setState(state);
        int action = policy->SelectAction();
        bool action_ok = environment->Act(action);
        if (!action_ok) {
            break;
        }
    }
    return statistics;
}
#endif
