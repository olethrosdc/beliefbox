/* -*- Mode: C++; -*- */
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "MountainCar.h"
#include "CartPole.h"
#include "Pendulum.h"
#include "MersenneTwister.h"
#include "DiscreteChain.h"
#include "EasyClock.h"
#include "BasisSet.h"
#include "RepresentativeStateValueIteration.h"
#include "ApproximateValueIteration.h"
#include "GaussianValueFunctionModel.h"
#include "KNNValueFunctionModel.h"
#include "debug.h"
#include <getopt.h>
#include <string>

struct Options {
    real gamma = 0.95;
    int grid_size = 4;
    real grid_scale = 1.0;
    int n_samples = 1;
    real threshold = 0;
    int n_iterations = 1000;
    int evaluation_grid_size = 32;
    std::string algorithm_name = "RSVI";
    std::string environment_name = "Pendulum";
	bool show_value_estimate = false;
	bool show_expected_utility = false;
	bool show_performance = false;
	int n_eval_samples = 1;
	int eval_horizon = 200; 
};

int main (int argc, char* argv[])
{
    MersenneTwisterRNG rng;
    Options options;
    {
        //options
        int c;
        while(1){
            int option_index = 0;
            static struct option long_options[] = {
                {"help", no_argument, 0, 0}, //0
                {"n_iterations",required_argument, 0, 0}, //1
                {"gamma", required_argument, 0, 0}, //2
                {"n_samples",required_argument, 0 ,0}, //3
                {"environment", required_argument, 0, 0}, //4
                {"grid", required_argument, 0, 0}, //5
                {"scale", required_argument, 0, 0}, //6
                {"algorithm", required_argument, 0, 0}, //7
                {"threshold", required_argument, 0, 0}, //8
                {"evaluation_grid", required_argument, 0, 0}, //9
				{"show_value_estimate", no_argument, 0, 0}, //10
				{"show_expected_utility", no_argument, 0, 0}, //11
				{"n_eval_samples", required_argument, 0, 0}, //12
				{"eval_horizon", required_argument, 0, 0}, //13
				{"show_performance", no_argument, 0, 0}, //14
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
                case 1: options.n_iterations = atoi(optarg); break;
                case 2: options.gamma = atof(optarg); break;
                case 3: options.n_samples = atoi(optarg); break;
                case 4: options.environment_name = optarg; break;
                case 5: options.grid_size = atoi(optarg); break;
                case 6: options.grid_scale = atof(optarg); break;
                case 7: options.algorithm_name = optarg; break;
                case 8: options.threshold = atof(optarg); break;
                case 9: options.evaluation_grid_size = atoi(optarg); break;
				case 10: options.show_value_estimate = true; break;
				case 11: options.show_expected_utility = true; break;
                case 12: options.n_eval_samples = atoi(optarg); break;
				case 13: options.eval_horizon = atoi(optarg); break;
				case 14: options.show_performance = true; break;
                default:
                    printf ("Usage options\n");
                    for (int i=0; long_options[i].name; i++) {
                        printf("--%s", long_options[i].name);
                        switch (long_options[i].has_arg) {
                        case no_argument: printf ("\n"); break;
                        case optional_argument: printf (" [arg]\n"); break;
                        case required_argument: printf (" arg\n"); break;
                        }
                    }
                    exit(0);
                    break;
                }
                break;
            default:
                printf ("Usage options\n");
                for (int i=0; long_options[i].name; i++) {
                    printf("--%s", long_options[i].name);
                    switch (long_options[i].has_arg) {
                    case no_argument: printf ("\n"); break;
                    case optional_argument: printf (" [arg]\n"); break;
                    case required_argument: printf (" arg\n"); break;
                    }
                }
                exit(-1);
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


    
    
    Environment<Vector, int>* environment_ptr  = NULL;
    if (options.environment_name == "MountainCar") {
        environment_ptr = new MountainCar;
    } else if (options.environment_name == "Pendulum") {
        environment_ptr = new Pendulum;
    } else if (options.environment_name == "CartPole") {
        environment_ptr = new CartPole;
    } else {
        std::cerr << "Unknown environment: " <<  options.environment_name << std::endl;
        std::cerr << "Choices: MountainCar, Pendulum, CartPole" << std::endl;
        exit(-1);
    }
    Environment<Vector, int>& environment = *environment_ptr;
                                          
    logmsg("Generating grid\n");
    EvenGrid grid(environment.StateLowerBound(),
                  environment.StateUpperBound(),
                  options.grid_size);

    logmsg("Creating kernel\n");
    // use an RBF basis for the kernel fromthe grid
    RBFBasisSet kernel(grid, options.grid_scale);

    logmsg("Selecting representative states\n");
    // create the set of representative states (identical to the grid)
    std::vector<Vector> states;
    std::vector<int> actions;
    for (int i=0; i<grid.getNIntervals(); ++i) {
        states.push_back(grid.getCenter(i));
    }
    // just use all the actions since it's a discrete space
    for (int i=0; i<environment.getNActions(); ++i) {
        actions.push_back(i);
    }

    ValueFunctionAlgorithm<Vector, int>* VFA = NULL;
	ValueFunctionModel<Vector, int>* VFM = NULL;
    if (options.algorithm_name == "RSVI") {
        logmsg("setting up RSVI\n");
        RepresentativeStateValueIteration<Vector, int, RBFBasisSet, Environment<Vector, int> >* rsvi = new RepresentativeStateValueIteration<Vector, int, RBFBasisSet, Environment<Vector, int> >(options.gamma, kernel, states, actions, environment, options.n_samples);
        rsvi->setThreshold(options.threshold);
        rsvi->setMaxIter(options.n_iterations);
		VFA = rsvi;
    } else if (options.algorithm_name == "AVI") {
        logmsg("setting up AVI\n");
		//VFM = new GaussianValueFunctionModel(environment.getNStates(), environment.getNActions());
		VFM = new KNNValueFunctionModel(environment.getNStates(), environment.getNActions(), 3);
		//VFM = new MultivariateNormalValueFunctionModel(environment.getNStates(), environment.getNActions());
        ApproximateValueIteration<Vector, int, ValueFunctionModel<Vector, int>, Environment<Vector, int> >* AVI = new ApproximateValueIteration<Vector, int, ValueFunctionModel<Vector, int>, Environment<Vector, int> >(options.gamma, *VFM, states, actions, environment, options.n_samples);
        AVI->setMaxIter(options.n_iterations);
		VFA = AVI;
    } else {
		
        std::cerr << "Unknown algorithm " << options.algorithm_name << std::endl;
        std::cerr << "Choices: RSVI, AVI"  << std::endl;
    }
    

	if (options.show_value_estimate) {
		logmsg("Calculating approximate value function\n");
		VFA->CalculateValues();
		std::string filename;
		filename.append(options.algorithm_name);
		filename.append("_"); filename.append(options.environment_name);
		filename.append("_"); filename.append(std::to_string(options.grid_size));
		filename.append("_"); filename.append(std::to_string(options.grid_scale));
		filename.append("_"); filename.append(std::to_string(options.n_samples));
		filename.append("_"); filename.append(std::to_string(options.n_iterations));
		filename.append(".values"); 
		FILE* outfile = fopen(filename.c_str(), "w");
		if (outfile) {
			EvenGrid evaluation_grid(environment.StateLowerBound(),
									 environment.StateUpperBound(),
									 options.evaluation_grid_size);
			for (int i=0; i<evaluation_grid.getNIntervals(); ++i) {
				Vector state = evaluation_grid.getCenter(i);
				state.printf(outfile);
				fprintf(outfile, " %f\n", VFA->getValue(state));
			}
			fclose(outfile);
		} else {
			Serror("Failed to write to file %s\n", filename.c_str());
		}
	}

	if (options.show_expected_utility) {
		logmsg("Estimating performance\n");
		std::string filename;
		filename.append(options.algorithm_name);
		filename.append("_"); filename.append(options.environment_name);
		filename.append("_"); filename.append(std::to_string(options.grid_size));
		filename.append("_"); filename.append(std::to_string(options.grid_scale));
		filename.append("_"); filename.append(std::to_string(options.n_samples));
		filename.append("_"); filename.append(std::to_string(options.n_iterations));
		filename.append(".utility"); 
		FILE* outfile = fopen(filename.c_str(), "w");
		if (outfile) {
			EvenGrid evaluation_grid(environment.StateLowerBound(),
									 environment.StateUpperBound(),
									 options.evaluation_grid_size);
			for (int i=0; i<evaluation_grid.getNIntervals(); ++i) {
				Vector start_state = evaluation_grid.getCenter(i);
				start_state.printf(outfile);
				real U = 0;
				for (int k=0; k<options.n_eval_samples; k++) {
					environment.Reset();
					environment.setState(start_state);
					real discount = 1;
					for (int t=0; t<options.eval_horizon; t++) {
						Vector state = environment.getState();
						real Q_max = -INF;
						real a_max = -1;

						for (int a=0; a<environment.getNActions(); ++a) {
							real Q_a = VFA->getValue(state, a);
							if (Q_a > Q_max) {
								Q_max = Q_a;
								a_max = a;
							}
						}
						bool action_ok = environment.Act(a_max);
						U += discount * environment.getReward();
						discount *= options.gamma;
						if (!action_ok) {
							break;
						}
					}
				}
				fprintf(outfile, " %f\n", U / (real) options.n_eval_samples);
			}
			fclose(outfile);
		} else {
			Serror("Failed to write to file %s\n", filename.c_str());
		}
	}

	if (options.show_performance) {
		logmsg("Estimating average performance\n");
		std::string filename;
		filename.append(options.algorithm_name);
		filename.append("_"); filename.append(options.environment_name);
		filename.append("_"); filename.append(std::to_string(options.grid_size));
		filename.append("_"); filename.append(std::to_string(options.grid_scale));
		filename.append("_"); filename.append(std::to_string(options.n_samples));
		filename.append("_"); filename.append(std::to_string(options.n_iterations));
		filename.append(".performance"); 
		FILE* outfile = fopen(filename.c_str(), "w");
		if (outfile) {
			EvenGrid evaluation_grid(environment.StateLowerBound(),
									 environment.StateUpperBound(),
									 options.evaluation_grid_size);
			VFA->setMaxIter(1);
			for (int iter=0; iter<options.n_iterations; iter++) {
				real U = 0;
				VFA->CalculateValues();
				for (int i=0; i<evaluation_grid.getNIntervals(); ++i) {
					Vector start_state = evaluation_grid.getCenter(i);
					for (int k=0; k<options.n_eval_samples; k++) {
						environment.Reset();
						environment.setState(start_state);
						real discount = 1;
						for (int t=0; t<options.eval_horizon; t++) {
							Vector state = environment.getState();
							real Q_max = -INF;
							real a_max = -1;
							
							for (int a=0; a<environment.getNActions(); ++a) {
								real Q_a = VFA->getValue(state, a);
								if (Q_a > Q_max) {
									Q_max = Q_a;
									a_max = a;
								}
							}
							bool action_ok = environment.Act(a_max);
							U += discount * environment.getReward();
							discount *= options.gamma;
							if (!action_ok) {
								break;
							}
						}
					}
				}
				U *= (real) (evaluation_grid.getNIntervals() * options.n_eval_samples);
				logmsg("U: %f\n", U);
				fprintf(outfile, " %f\n", U);
				
			}
			fclose(outfile);
		} else {
			Serror("Failed to write to file %s\n", filename.c_str());
		}
	}
	
	printf("\nDone\n");
	delete environment_ptr;
	delete VFM;
	delete VFA;
	return 0.0;
}


#endif
