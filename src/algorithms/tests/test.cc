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
#include "Demonstrations.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "RandomPolicy.h"
#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "KNNClassifier.h"
#include "ClassifierMixture.h"
#include "ClassifierPolicy.h"
#include "EasyClock.h"
#include <cstring>
#include <getopt.h>

template <class X, class A>
class TotalRewardStatistic
{
public:
    real getAverageTotalReward(Demonstrations<X, A>& data)
    {
        real r_d = 0;
        for (uint i=0; i<data.total_rewards.size(); ++i) {
            r_d += data.total_rewards[i];
        }
        r_d /= (real) data.total_rewards.size();
        return r_d;
    }

    real distance(Demonstrations<X, A>& data,
                  Demonstrations<X, A>& sample)
    {
        real r_d = getAverageTotalReward(data);
        real r_s = getAverageTotalReward(sample);
		//printf ("%f %f ", r_d, r_s);
        return fabs(r_d - r_s);
    }

};


/** This template has 5 types.

 G: The generator class, which must provide
 - M G::generate()

 F: The statistic class, which must provide

 M: The model class, which must provide

 P: The policy class which must provide

 X: The observation class
 
 A: the action

*/
template <class G, class F, class M, class P, typename X, typename A>
class ABCRL
{
public:
	G generator;
	F statistic;
	ABCRL()
	{
	}
	M GetSample(Demonstrations<X,A>& data,
                P& policy,
                real discounting,
                real epsilon,
                int max_iter)
	{
		//real min_epsilon = INF;
		for (int iter=0; iter<max_iter; ++iter) {
            M model = generator.generate();
            Demonstrations<X, A> sample;
            for (int i=0; i<100; ++i) {
                sample.Simulate(model, policy, discounting, -1);
            }
            real error = statistic.distance(data, sample);
            //printf("%f # error\n", error);
            if (error <= epsilon) {
				printf ("%f ", error); 
				model.Show();
                //return model;
            }
		}
        Swarning("No model generated\n");
        return generator.generate();
	}
	
};

static const char* const help_text = "Usage: test [options]\n\
\nOptions:\n\
    --environment: {MountainCar, Pendulum}\n\
    --gamma:       reward discounting in [0,1]\n\
\n";

int main(int argc, char* argv[])
{
    
	// Create a new environment
	Environment<Vector, int>* environment;

	MersenneTwisterRNG rng;

	struct Options
	{
		real gamma;
        real epsilon;
		char* environment_name;
	};

	Options options;
	options.gamma = 0.99;
	options.epsilon = 10.0;
	options.environment_name = NULL;

	{
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"gamma", required_argument, 0, 0}, //0
                {"environment", required_argument, 0, 0}, //1
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
                case 0: options.gamma = atof(optarg); break;
                case 1: options.environment_name = optarg; break;
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

    
    if (!options.environment_name) {
        fprintf(stderr, "Must specify environment\n");
        exit(-1);
    }
    if (!strcmp(options.environment_name, "MountainCar")) {
        environment = new MountainCar();
    } else {
        fprintf(stderr, "Invalid environment name %s\n", options.environment_name);
        exit(-1);
    }

	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;
	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), &rng);
    //template <class G, class F, class M, class P, typename X, typename A>
    MountainCar mountain_car(true);
    Demonstrations<Vector, int> data;

	mountain_car.Show();
    for (int i=0; i<100; ++i) {
        data.Simulate(mountain_car, *policy, options.gamma, -1);
    }

    ABCRL<MountainCarGenerator, TotalRewardStatistic<Vector, int>, MountainCar, AbstractPolicy<Vector, int>, Vector, int> abcrl;


    MountainCar sample = abcrl.GetSample(data,
                                         *policy,
                                         options.gamma,
                                         options.epsilon,
                                         1000);
	return 0;
}
#endif
