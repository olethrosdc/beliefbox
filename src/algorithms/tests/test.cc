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
#include "LSTDQ.h"
#include "LSPI.h"
#include "Grid.h"
#include "BasisSet.h"

#include "RSAPI.h"
#include "Rollout.h"
#include "Demonstrations.h"

#include "MountainCar.h"
#include "Pendulum.h"
#include "PuddleWorld.h"

#include "RandomPolicy.h"

#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"

#include "KNNClassifier.h"
#include "ClassifierMixture.h"
#include "ClassifierPolicy.h"
#include "EasyClock.h"

#include <cstring>
#include <getopt.h>
#include <vector>

/** Options */
struct Options
{
    real gamma; ///< discount factor
    real epsilon; ///< stopping threshold
    char* environment_name; ///< environment name
    RandomNumberGenerator& rng; ///< random number generator
    int n_trajectories; ///< number of trajectories to take for each sample
    int n_samples; ///< number of sampled environments
    int n_training; ///< number of training trajectories
    int n_testing; ///< number of testing trajectories
    int grid; ///< grid size for features
    real scale; ///< scale of RBF basis
    real accuracy; ///< value/policy iteration accuracy
    int lspi_iterations; ///< number of lspi iterations
    int n_evaluations; ///< number of evaluations
    bool reuse_training_data; ///< reuse training data in lspi
};


/** The total reward statistic
    
    X - observation
    A - action
*/
template <class X, class A>
class TotalRewardStatistic
{
protected:
    real delta;
    real C;
public:
    TotalRewardStatistic(real delta_ = 0.5)
        : delta(delta_), C(0.5*log(2.0/delta))
    {
    }
    
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
		//printf ("%f %f (r)", r_d, r_s);
        
        return fabs(r_d - r_s) + sqrt(C * (1.0 / (real) data.size() + 1.0 / (real) sample.size()));
    }
};

/** get LSPI policy

    If the environment is given, then always rollout from the environment.
    If the data is given, and options.use_training_data is true,
    then add these rollouts to the original data.
 */
AbstractPolicy<Vector, int>*
getLSPIPolicy(Environment<Vector, int>* training_environment,
              AbstractPolicy<Vector, int>& policy,
              Rollout<Vector, int, AbstractPolicy<Vector, int> >* data,
              RBFBasisSet& RBFs,
              Options& options)
{
    Environment<Vector, int>* environment;
    if (training_environment) {
        environment = training_environment;
    } else {
        environment = data->environment;
    }
    int state_dimension = environment->getNStates();
    int n_actions = environment->getNActions();
    Vector S_L = environment->StateLowerBound();
    Vector S_U = environment->StateUpperBound();

    Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
    if (data && options.reuse_training_data) {
        // copy the old data in the new rollout!
        rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> > (*data);
    } else {
        rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> > (urandom(S_L, S_U), &policy, environment, options.gamma, true);
    }
    if (training_environment) {
        rollout->UniformSampling(options.n_training, -1);
    }
    
    logmsg("Total number of collected rollout samples -> %d\n", rollout->getNSamples());
    EvenGrid Discretisation(S_L, S_U, options.grid);
    LSPI lspi(options.gamma, options.accuracy, state_dimension, n_actions, options.lspi_iterations, &RBFs, rollout);
    lspi.PolicyIteration();
    AbstractPolicy<Vector, int> * lspi_policy = new FixedContinuousPolicy(lspi.ReturnPolicy());
    delete rollout;
    return lspi_policy;
}

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
	void GetSample(Demonstrations<X,A>& data,
                   P& policy,
                   real discounting,
                   real epsilon,
                   int n_trajectories,
                   int n_samples,
                   std::vector<M>& samples,
                   std::vector<real>& values)
	{
		//real min_epsilon = INF;
        int n_models = 0;
		for (int iter=0; iter<n_samples || n_models == 0; ++iter) {
            M model = generator.Generate();
            Demonstrations<X, A> sample;
            for (int i=0; i<n_trajectories; ++i) {
                sample.Simulate(model, policy, discounting, -1);
            }
            real error = statistic.distance(data, sample);
            if (error <= epsilon) {
                n_models++;
                samples.push_back(model);
                values.push_back(statistic.getAverageTotalReward(sample));
				printf ("# %d (e) %f %f\n", n_models, error, epsilon); 
				model.Show();
            }
			if(iter > n_samples) {
				iter = 0;
				epsilon *= 2.0;
				printf("# e -> %f\n", epsilon);
				fflush(stdout);
			}
		}
        if (n_models == 0) {
			Swarning("No model generated: %d\n", (int) samples.size());
		}
	}
	
};


template <class X, class A, class M, class P>
real EvaluatePolicy(M& environment, P& policy, real gamma, int n_testing)
{
	int n_samples = 0;
	real discounted_reward = 0;
    real average_steps = 0;
    for (int i=0; i<n_testing; ++i) {
        Demonstrations<Vector, int> data;
        data.Simulate(environment, policy, gamma, -1);
		for (uint t=0; t<data.size(); ++t) {
            average_steps += data.length(t);
			discounted_reward += data.discounted_rewards[t];
			n_samples++;
		}
    }
    logmsg("Average steps: %f\n", average_steps / (real) n_samples);
	return discounted_reward / (real) n_samples;
}

/** Run a test
 */
template <class G, class M>
void RunTest(Options& options)
{

    G generator;
    M environment = generator.Generate();



    int state_dimension = environment.getNStates();
    //int n_actions = environment.getNActions();
    Vector S_L = environment.StateLowerBound();
    Vector S_U = environment.StateUpperBound();

    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
	 
    EvenGrid Discretisation(S_L, S_U, options.grid);
    RBFBasisSet RBFs(Discretisation, options.scale);




	// Start with a random policy!
	RandomPolicy random_policy(environment.getNActions(), &options.rng);

	// Placeholder for the policy
	AbstractPolicy<Vector, int>& policy = random_policy;

    //template <class G, class F, class M, class P, typename X, typename A>
    Demonstrations<Vector, int> training_data;
	logmsg("Training\n"); environment.Show();
    for (int i=0; i<options.n_training; ++i) {
        training_data.Simulate(environment, policy, options.gamma, -1);
    }
    
    Rollout<Vector, int, AbstractPolicy<Vector, int> > training_rollouts(urandom(S_L, S_U), &policy, &environment, options.gamma, true);
    training_rollouts.StartingDistributionSampling(options.n_training, -1);


    ABCRL<G, TotalRewardStatistic<Vector, int>, M, AbstractPolicy<Vector, int>, Vector, int> abcrl;


    std::vector<M> samples;
    std::vector<real> values;
    abcrl.GetSample(training_data,
                    policy,
                    options.gamma,
                    options.epsilon,
                    options.n_trajectories,
                    options.n_samples,
                    samples,
                    values);

	real V_initial = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, policy, options.gamma, options.n_testing);
	Vector V((uint) samples.size());
	Vector hV((uint) samples.size());
	Vector V_LSPI((uint) samples.size());
	Vector hV_LSPI((uint) samples.size());
	for (int i=0; i<V.Size(); ++i) {
		V(i) = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(samples[i], policy, options.gamma, options.n_testing);
		hV(i) = values[i];

        AbstractPolicy<Vector, int>* lspi_policy
            =  getLSPIPolicy(&samples[i],
                             policy,
                             &training_rollouts,
                             RBFs,
                             options);
        
        hV_LSPI(i)= EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(samples[i], *lspi_policy, options.gamma, options.n_testing);

        V_LSPI(i) = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *lspi_policy, options.gamma, options.n_testing);
		printf ("%f %f %f %f# sampled value\n", hV(i), V(i), hV_LSPI(i), V_LSPI(i));
        fflush(stdout);
	}

    AbstractPolicy<Vector, int>* oracle_lspi_policy
        =  getLSPIPolicy(NULL,
                         policy,
                         &training_rollouts,
                         RBFs,
                         options);
    real V_lspi_oracle = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *oracle_lspi_policy, options.gamma, options.n_testing);
    real V_LSTD = EvaluateLSTD(environment, policy, training_data, options);    
    printf ("%f %f %f %f %f %f %f # hV V V_oracle hV_LSPI V_LSPI V_LSPI_oracle V_LSTD\n",
			hV.Sum()/(real) hV.Size(),
			V.Sum()/(real) V.Size(),
            V_initial,
            hV_LSPI.Sum() / (real) hV_LSPI.Size(),
            V_LSPI.Sum() / (real) V_LSPI.Size(),
            V_lspi_oracle,
            V_LSTD);
}

/** Run an online test
 */
template <class G, class M>
void RunOnlineTest(Options& options)
{
    
    // First, generate the environment
    G generator;
    M environment = generator.Generate();

    int state_dimension = environment.getNStates();
    //int n_actions = environment.getNActions();
    Vector S_L = environment.StateLowerBound();
    Vector S_U = environment.StateUpperBound();

    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
	 
    // Discretise it according to the grid scale
    EvenGrid Discretisation(S_L, S_U, options.grid);
    RBFBasisSet RBFs(Discretisation, options.scale);

	// Start with a random policy!
	RandomPolicy random_policy(environment.getNActions(), &options.rng);

    // Policy placeholder
	AbstractPolicy<Vector, int>& policy = random_policy;

    //template <class G, class F, class M, class P, typename X, typename A>
    Demonstrations<Vector, int> training_data;
	logmsg("Training\n"); environment.Show();
    for (int i=0; i<options.n_training; ++i) {
        training_data.Simulate(environment, policy, options.gamma, -1);
    }
    
    Rollout<Vector, int, AbstractPolicy<Vector, int> > training_rollouts(urandom(S_L, S_U), &policy, &environment, options.gamma, true);
    training_rollouts.StartingDistributionSampling(options.n_training, -1);


    ABCRL<G, TotalRewardStatistic<Vector, int>, M, AbstractPolicy<Vector, int>, Vector, int> abcrl;


    std::vector<M> samples;
    std::vector<real> values;
    abcrl.GetSample(training_data,
                    policy,
                    options.gamma,
                    options.epsilon,
                    options.n_trajectories,
                    options.n_samples,
                    samples,
                    values);

	real V_initial = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, policy, options.gamma, options.n_testing);
	Vector V((uint) samples.size());
	Vector hV((uint) samples.size());
	Vector V_LSPI((uint) samples.size());
	Vector hV_LSPI((uint) samples.size());
	for (int i=0; i<V.Size(); ++i) {
		V(i) = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(samples[i], policy, options.gamma, options.n_testing);
		hV(i) = values[i];

        AbstractPolicy<Vector, int>* lspi_policy
            =  getLSPIPolicy(&samples[i],
                             policy,
                             &training_rollouts,
                             RBFs,
                             options);
        
        hV_LSPI(i)= EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(samples[i], *lspi_policy, options.gamma, options.n_testing);

        V_LSPI(i) = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *lspi_policy, options.gamma, options.n_testing);
		printf ("%f %f %f %f# sampled value\n", hV(i), V(i), hV_LSPI(i), V_LSPI(i));
        fflush(stdout);
	}

    AbstractPolicy<Vector, int>* oracle_lspi_policy
        =  getLSPIPolicy(NULL,
                         policy,
                         &training_rollouts,
                         RBFs,
                         options);
    real V_lspi_oracle = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *oracle_lspi_policy, options.gamma, options.n_testing);
    real V_LSTD = EvaluateLSTD(environment, policy, training_data, options);    
    printf ("%f %f %f %f %f %f %f # hV V V_oracle hV_LSPI V_LSPI V_LSPI_oracle V_LSTD\n",
			hV.Sum()/(real) hV.Size(),
			V.Sum()/(real) V.Size(),
            V_initial,
            hV_LSPI.Sum() / (real) hV_LSPI.Size(),
            V_LSPI.Sum() / (real) V_LSPI.Size(),
            V_lspi_oracle,
            V_LSTD);
}


real EvaluateLSTD(Environment<Vector, int>& environment,
				  AbstractPolicy<Vector, int>& policy,
				  Demonstrations<Vector, int>& data,
                  Options& options)
{
    int state_dimension = environment.getNStates();
    int n_actions = environment.getNActions();
    Vector S_L = environment.StateLowerBound();
    Vector S_U = environment.StateUpperBound();

    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
	 
    EvenGrid Discretisation(S_L, S_U, options.grid);
    RBFBasisSet RBFs(Discretisation, options.scale);
	LSTDQ lstdq(options.gamma,
				state_dimension,
				n_actions,
				RBFs,
				data);
	
	//lstdq.Calculate_Opt();
	lstdq.Calculate();
	
	real V = 0;
	for (int i=0; i<options.n_testing; ++i) {
		environment.Reset();
		Vector state = environment.getState();
		policy.Reset();
		policy.setState(state);
		real Vi = lstdq.getValue(state, rand()%n_actions);
		//printf ("%f ", Vi);
		V += Vi;
	}
	//printf("# LSTD values\n");
	V /= (real) options.n_testing;
	return V;
}



static const char* const help_text = "Usage: test [options]\n\
\nOptions:\n\
    --environment:           {MountainCar, Pendulum, Puddleworld}\n\
    --discount:              reward discounting in [0,1]\n\
    --threshold:             statistic threshold\n\
    --n_trajectories:        number of trajectories per sample\n\
    --n_samples:             number of sampled models\n\
    --n_training:            number of training trajectories\n\
    --n_testing:             number of test trajectories\n\
    --seed:                  seed all the RNGs with this\n\
    --grid:                  number of grid intervals for LSTD\n\
    --scale:                 RBF scale for LSTD\n\
    --n_evaluations:         number of evaluations\n\
    --reuse_training_data:   reuse the training data in LSPI samples\n\
    --online:                do the online test\n\
\n";

int main(int argc, char* argv[])
{
	int seed = time(NULL);
    bool online_test = false;

	MersenneTwisterRNG rng;
	Options options = 
        {0.99, // discount
         10.0, // threshold
         NULL, // envionment
         rng, // rng
         1000, // trajectories per sample
         128, // samples
         10, // training trajectories
         10000, // test trajectories
         5, // grid size
         0.5, // rbf scale
         1e-6, // accuracy
         100, // lspi iterations
         1000 // number of evaluations
        };

	{
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"discount", required_argument, 0, 0}, //0
                {"environment", required_argument, 0, 0}, //1
                {"threshold", required_argument, 0, 0}, //2
                {"n_trajectories", required_argument, 0, 0}, //3
                {"n_samples", required_argument, 0, 0}, //4
                {"n_training", required_argument, 0, 0}, //5
                {"n_testing", required_argument, 0, 0}, //6
                {"grid", required_argument, 0, 0}, //7
                {"scale", required_argument, 0, 0}, //8
				{"seed", required_argument, 0, 0}, //9
				{"accuracy", required_argument, 0, 0}, //10
				{"lspi_iterations", required_argument, 0, 0}, //11
				{"n_evaluations", required_argument, 0, 0}, //12
                {"reuse_training_data", no_argument, 0, 0}, //13
                {"online", no_argument, 0, 0}, //14
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
                case 2: options.epsilon = atof(optarg); break;
                case 3: options.n_trajectories = atoi(optarg); break;
                case 4: options.n_samples = atoi(optarg); break;
                case 5: options.n_training = atoi(optarg); break;
                case 6: options.n_testing = atoi(optarg); break;
                case 7: options.grid = atoi(optarg); break;
                case 8: options.scale = atof(optarg); break;
				case 9: seed = atoi(optarg); break;
                case 10: options.accuracy = atof(optarg); break;
                case 11: options.lspi_iterations = atoi(optarg); break;
                case 12: options.n_evaluations = atoi(optarg); break;
                case 13: options.reuse_training_data = true; break;
                case 14: online_test = true; break;
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

	srand(seed);
    srand48(seed);
    rng.manualSeed(seed);
	setRandomSeed(seed);

    if (!options.environment_name) {
        Serror("Must specify environment\n");
        exit(-1);
    }
    if (options.gamma < 0 || options.gamma > 1) {
        Serror("gamma must be in [0,1]\n");
        exit(-1);
    }

    if (options.n_samples < 1) {
        Serror("n_samples must be >= 1\n");
        exit(-1);
    }

    if (options.n_trajectories < 1) {
        Serror("n_trajectories must be >= 1\n");
        exit(-1);
    }

    if (options.n_training < 0) {
        Serror("n_training must be >= 0\n");
        exit(-1);
    }

    if (options.n_testing < 1) {
        Serror("n_training must be >= 1\n");
        exit(-1);
    }

    logmsg("Starting environment %s\n", options.environment_name);

    if (online_test) {
        for (int i=0; i<options.n_evaluations; ++i) {
            if (!strcmp(options.environment_name, "MountainCar")) {
            logmsg("Testing mountain car\n");
            RunOnlineTest<MountainCarGenerator, MountainCar>(options);
            } else if (!strcmp(options.environment_name, "Pendulum")) {
                logmsg("Testing pendulum\n");
                RunOnlineTest<PendulumGenerator, Pendulum>(options);
            } else if (!strcmp(options.environment_name, "PuddleWorld")) {
                logmsg("Testing puddle world\n");
                RunOnlineTest<PuddleWorldGenerator, PuddleWorld>(options);
            } else {
                fprintf(stderr, "Invalid environment name %s\n", options.environment_name);
                exit(-1);
            }
        }
    } else {
        for (int i=0; i<options.n_evaluations; ++i) {
            if (!strcmp(options.environment_name, "MountainCar")) {
            logmsg("Testing mountain car\n");
            RunTest<MountainCarGenerator, MountainCar>(options);
            } else if (!strcmp(options.environment_name, "Pendulum")) {
                logmsg("Testing pendulum\n");
                RunTest<PendulumGenerator, Pendulum>(options);
            } else if (!strcmp(options.environment_name, "PuddleWorld")) {
                logmsg("Testing puddle world\n");
                RunTest<PuddleWorldGenerator, PuddleWorld>(options);
            } else {
                fprintf(stderr, "Invalid environment name %s\n", options.environment_name);
                exit(-1);
            }
        }
    }
	return 0;
}
#endif
