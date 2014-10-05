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


#include "UCTMC.h"
#include "UCTMCL.h"


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
#include "Bike.h"
#include "Acrobot.h"
#include "CartPole.h"

#include "RandomPolicy.h"

#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFile.h"

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
  int n_model_rollouts; ///< number of rollouts to take from the model
  int grid; ///< grid size for features
  real scale; ///< scale of RBF basis
  real accuracy; ///< value/policy iteration accuracy
  int lspi_iterations; ///< number of lspi iterations
  int n_evaluations; ///< number of evaluations
  bool reuse_training_data; ///< reuse training data in lspi
  int n_episodes; ///< number of episodes
  bool sampling; ///< use sampling online
  real delta; ///< error probability bound
  real R_max; ///< maximum reward
  real c_uct; ///< uct factor
  int n_rollouts_uct; ///< number of rollouts for uct
  real lambda_uct; ///<  lambda mixing for uct
  real stepsize_uct; ///< step size
  int depth_uct; ///< maximum tree depth
  int grid_uct; ///< grid utc
  Options(RandomNumberGenerator& rng_) :
    gamma(0.99),
    epsilon(1.0),
    environment_name(NULL),
    rng(rng_), n_trajectories(1000),
    n_samples(128),
    n_training(10),
    n_testing(1000),
    n_model_rollouts(2000),
    grid(5),
    scale(0.5),
    accuracy(1e-6),
    lspi_iterations(25),
    n_evaluations(100),
    reuse_training_data(false),
    sampling(false),
    delta(0.5),
    R_max(1.0),
    c_uct(1000),
    n_rollouts_uct(1000),
    lambda_uct(1),
    stepsize_uct(0.001),
    depth_uct(1000),
    grid_uct(20)
  {
  }
     
  void ShowOptions()
  {
    logmsg("------------------------\n");
    logmsg("Options\n");
    logmsg("=======\n");
    logmsg("Gamma: %f\n", gamma);
    logmsg("Epsilon: %f\n", epsilon);
    logmsg("Trajectories: %d\n", n_trajectories);
    logmsg("n_samples: %d\n", n_samples);
    logmsg("n_training: %d\n",  n_training); 
    logmsg("n_testing: %d\n", n_testing);
    logmsg("n_model_rollouts: %d\n", n_model_rollouts);
    logmsg("grid %d\n", grid);
    logmsg("scale %f\n", scale);
    logmsg("accuracy %f\n", accuracy);
    logmsg("lspi iterations %d\n", lspi_iterations);
    logmsg("n_evaluations %d\n", n_evaluations);
    logmsg("reuse_training_data %d\n", reuse_training_data);
    logmsg("n_episodes %d\n", n_episodes);
    logmsg("sampling %d\n", sampling);
    logmsg("delta %f\n", delta);
    logmsg("R_max %f\n", R_max);
    logmsg("c_uct %f\n", c_uct);
    logmsg("n_rollouts_uct %d\n", n_rollouts_uct);
    logmsg("lambda_uct %f\n", lambda_uct);
    logmsg("stepsize_uct %f\n", stepsize_uct);
    logmsg("depth_uct %d\n", depth_uct); 
    logmsg("grid_uct %d\n", grid_uct); 
    logmsg("------------------------\n");
  }
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
  TotalRewardStatistic(real delta_ = 1.0)
    : delta(delta_), C(0.5*log(1.0/delta))
  {
  }
  void setBound(real delta_, real gamma, real R_max)
  {
    delta = delta_;
    C = 0.5 * R_max * R_max * log(1.0/delta) / (1.0 - gamma);
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
    real n_d = (real) data.size();
    real n_s = (real) sample.size();
#if 1
    real bound = sqrt(C * (1.0 / n_d + 1.0 / n_s));
    real distance = fabs(r_d - r_s) - bound;
    //logmsg ("statistic: %f %f (r) %f | %f => %f \n", r_d, r_s, bound, fabs(r_d - r_s), distance);
#else
    real gap = r_d - r_s;
    real distance = fabs(gap); //+ exp( - gap*gap / (1.0/n_d + 1.0/n_s));
    //logmsg ("statistic: %f %f \n", gap, distance);
#endif
    return distance; 
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
  if (data && (!training_environment || options.reuse_training_data)) {
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
  void GetSample(Demonstrations<X,A>& data, ///< input data
		 P& policy, ///< input policy
		 real discounting, ///< discounting
		 real epsilon, ///< acceptance margin
		 real delta, ///< error probability
		 real R_max, ///< reward bound
		 int n_trajectories, ///< number of trajectories to draw
		 int n_samples, ///< number of samples to generate
		 std::vector<M>& samples, ///< model samples
		 std::vector<real>& values)
  {
    //real min_epsilon = INF;
    int n_models = 0;
    statistic.setBound(delta, discounting, R_max);
    for (int iter=0; n_models == 0; ++iter) {
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
	logmsg ("models: %d (e) %f %f\n", n_models, error, epsilon); 
	model.Show();
      }
      if(iter > n_samples) {
	iter = 0;
	epsilon *= 2.0;
	logmsg("# increasing epsilon -> %f\n", epsilon);
	fflush(stdout);
      }
    }
    if (n_models == 0) {
      Swarning("No model generated: %d\n", (int) samples.size());
    }
  }
  void GetSample(std::vector<Demonstrations<X,A> >& data_list,
		 std::vector<P*>& policy_list,
		 real discounting,
		 real epsilon,
		 real delta,
		 real R_max,
		 int n_trajectories,
		 int n_samples,
		 std::vector<M>& samples,
		 std::vector<real>& values)
  {
    assert(policy_list.size() == data_list.size());

    int list_size = data_list.size();
    //real min_epsilon = INF;
    int n_models = 0;
    statistic.setBound(delta, discounting, R_max);
    for (int iter=0; n_models == 0; ++iter) {
      M model = generator.Generate();
      real error = 0.0;
      for (int k=0; k<list_size; ++k) {
	Demonstrations<X, A> sample;
	for (int i=0; i<n_trajectories; ++i) {
	  sample.Simulate(model, *policy_list[k], discounting, -1);
	}
	error += statistic.distance(data_list[k], sample);
      }
      error /= (real) list_size;
      //logmsg("total error: %f \n", error);
      if (error <= epsilon) {
	n_models++;
	samples.push_back(model);
	model.Show();
	logmsg ("models: %d (e) %f %f\n", n_models, error, epsilon); 
      }
      if(iter > n_samples) {
	iter = 0;
	epsilon *= 2.0;
	logmsg("epsilon increased -> %f\n", epsilon);
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
  int horizon =  (int) ceil(10.0/(1.0 - gamma));
  logmsg("Evaluating policy with horizon %d\n", horizon);
  Demonstrations<Vector, int> data;
  for (int i=0; i<n_testing; ++i) {
    data.Simulate(environment, policy, gamma, horizon);
  }
  for (uint t=0; t<data.size(); ++t) {
    average_steps += data.length(t);
    discounted_reward += data.discounted_rewards[t];
    //discounted_reward += data.total_rewards[t];
    n_samples++;
  }
  logmsg("Average steps: %f\n", average_steps / (real) n_samples);
  return discounted_reward / (real) n_samples;
}

/** Run a test
 */
template <class G, class M>
void RunTest(Options& options)
{

  /// Initialise the generator and generate a FIXED environment
  G generator;
  M environment = generator.Generate(false);

  int state_dimension = environment.getNStates();
  //int n_actions = environment.getNActions();
  Vector S_L = environment.StateLowerBound();
  Vector S_U = environment.StateUpperBound();

  logmsg("State dimension: %d\n", state_dimension);
  logmsg("S_L: "); S_L.print(stdout);
  logmsg("S_U: "); S_U.print(stdout);
	 
  EvenGrid Discretisation(S_L, S_U, options.grid);
  RBFBasisSet RBFs(Discretisation, options.scale);

  // Start with a random policy!
  RandomPolicy random_policy(environment.getNActions(), &options.rng);

  // Placeholder for the policy
  AbstractPolicy<Vector, int>& policy = random_policy;

  // generate demonstrations
  Demonstrations<Vector, int> training_data;
  logmsg("Training\n"); environment.Show();
  for (int i=0; i<options.n_training; ++i) {
    training_data.Simulate(environment, policy, options.gamma, -1);
  }
    
  // generate training rollouts 
  Rollout<Vector, int, AbstractPolicy<Vector, int> > training_rollouts(urandom(S_L, S_U), &policy, &environment, options.gamma, true);
  training_rollouts.StartingDistributionSampling(options.n_training, -1);

  // start ABC-RL
  ABCRL<G, TotalRewardStatistic<Vector, int>, M, AbstractPolicy<Vector, int>, Vector, int> abcrl;
    
  // start UCT
  EvenGrid discretize(environment.StateLowerBound(),environment.StateUpperBound(),options.grid_uct);	
  UCTMC<Vector, int> mcts(options.gamma, options.c_uct, &environment, &options.rng, discretize, options.stepsize_uct, options.lambda_uct, options.depth_uct, options.n_rollouts_uct);

  logmsg("ABC Generating samples\n");
  std::vector<M> samples;
  std::vector<real> values;
  double abc_time = 0;
  double abc_start = GetCPU();
  abcrl.GetSample(training_data,
		  policy,
		  options.gamma,
		  options.epsilon,
		  options.delta,
		  options.R_max,
		  options.n_trajectories,
		  options.n_samples,
		  samples,
		  values);
  abc_time += GetCPU() - abc_start;

  logmsg("Evaluating initial policy on %d samples\n", (int) samples.size());
  real V_initial = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, policy, options.gamma, options.n_testing);
    
  real V_ABC_UCT = 0;
  double abc_uct_time = 0;

  Vector V_LSPI((uint) samples.size());
  Options estimation_options = options;
  estimation_options.n_training = options.n_model_rollouts;
  logmsg("Running estimation policy with %d simulated trajectories\n", estimation_options.n_training);
  for (int i=0; i<V_LSPI.Size(); ++i) {
    logmsg("Estimating policy from simulation\n");
    double abc_start = GetCPU();
    AbstractPolicy<Vector, int>* lspi_policy
      =  getLSPIPolicy(&samples[i],
		       policy,
		       NULL,
		       RBFs,
		       estimation_options);
    abc_time += GetCPU() - abc_start;
    logmsg("Evaluating sampled policy\n");

    V_LSPI(i) = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *lspi_policy, options.gamma, options.n_testing);
    printf ("%f# sampled value\n", V_LSPI(i));
    fflush(stdout);
  }


  logmsg("Evaluating UCTMC policy\n");
  if (samples.size() > 0) {
    bool running;
    real reward;

    ContinuousStateEnvironment* sampled_environment = &samples[0];
    mcts.setEnvironment(sampled_environment);

    int state_dimension = sampled_environment->getNStates();
    Vector S_L = sampled_environment->StateLowerBound();
    Vector S_U = sampled_environment->StateUpperBound();
		
    logmsg("State dimension: %d\n", state_dimension);
    logmsg("S_L: "); S_L.print(stdout);
    logmsg("S_U: "); S_U.print(stdout);
		
    int horizon =  (int) ceil(10.0/(1.0 - options.gamma));
    V_ABC_UCT = 0;
    abc_uct_time = 0;
    double abc_uct_start = GetCPU();
    /// Run the UCT planner on the real thing
    for(int episode = 0; episode<options.n_testing; ++episode) {
      int step               = 0;
      real discounted_reward = 0;
      real total_reward      = 0;
      environment.Reset();
			
      Vector state = environment.getState();
      //state.print(stdout);
      int action =  mcts.PlanPolicy(state);
      do {
				
	running = environment.Act(action);
	reward = environment.getReward();
				
	total_reward += reward;
	discounted_reward += pow(options.gamma,step)*reward;
				
	state = environment.getState();
	//state.print(stdout);
	action = mcts.PlanPolicy(state);
	step++;
      } while(running && step <  horizon);
      //logmsg("Sampled Episode = %d: Steps = %d, Total Reward = %f, Discounted Reward = %f\n",episode, step, total_reward, discounted_reward);
      V_ABC_UCT += discounted_reward;
    }
    V_ABC_UCT /= (real) options.n_testing;
			
    abc_uct_time += GetCPU() - abc_uct_start;
  }
	
	
  double lspi_time = 0;
  double lspi_start = GetCPU();
  AbstractPolicy<Vector, int>* oracle_lspi_policy
    =  getLSPIPolicy(NULL,
		     policy,
		     &training_rollouts,
		     RBFs,
		     options);
  lspi_time += GetCPU() - lspi_start;
  logmsg("Evaluating LSPI policy\n");

  real V_lspi_oracle = EvaluatePolicy<Vector, int, M, AbstractPolicy<Vector, int> >(environment, *oracle_lspi_policy, options.gamma, options.n_testing);
    
  printf ("%f %f %f %f %f %f %f # V V_ABC_LSPI V_LSPI V_ABC_UCT T_ABC T_LSPI T_ABC_UCT\n",
	  V_initial,
	  V_LSPI.Sum() / (real) V_LSPI.Size(),
	  V_lspi_oracle,
	  V_ABC_UCT,
	  abc_time,
	  lspi_time,
	  abc_uct_time);
}

/** Run an online test
 */
template <class G, class M>
void RunOnlineTest(Options& options)
{
  int n_episodes_per_iteration = 1;
  int rollout_horizon = 1000;

  // First, generate the environment
  G generator;
  M environment = generator.Generate(false);
  logmsg("generated environment: ");
  environment.Show();
    
  int state_dimension = environment.getNStates();
  //int n_actions = environment.getNActions();
  Vector S_L = environment.StateLowerBound();
  Vector S_U = environment.StateUpperBound();

  logmsg("Running online test\n");
  logmsg("State dimension: %d\n", state_dimension);
  logmsg("S_L: "); S_L.print(stdout);
  logmsg("S_U: "); S_U.print(stdout);
	 
  // Discretise it according to the grid scale
  EvenGrid Discretisation(S_L, S_U, options.grid);
  RBFBasisSet RBFs(Discretisation, options.scale);


  std::vector<AbstractPolicy<Vector, int>* > policy_list;   
  std::vector<Demonstrations<Vector, int> > training_data_list;
    
  // Add the random policy to the list.
  policy_list.push_back(new RandomPolicy (environment.getNActions(), &options.rng));

  Rollout<Vector, int, AbstractPolicy<Vector, int> > training_rollouts(urandom(S_L, S_U), policy_list.back(), &environment, options.gamma, true);
    
  real epsilon_greedy = 1.0;
  logmsg("Running for %d episodes\n", options.n_training);
  for (int episode=0; episode < options.n_training; ++episode) {
    logmsg("Episode %d\n", episode); 
    fflush(stdout);
    AbstractPolicy<Vector, int>& policy = *(policy_list.back());
    if (options.sampling) {
      // for sampling
      training_data_list.push_back(Demonstrations<Vector, int> ());
      Demonstrations<Vector, int>& training_data = training_data_list.back();
      for (int i=0; i<n_episodes_per_iteration; ++i) {
	training_data.Simulate(environment, policy, options.gamma, rollout_horizon);
      }
      printf("%f %f %d # total, discounted reward, steps ABC\n",
	     training_data.total_rewards.back(),
	     training_data.discounted_rewards.back(),
	     training_data.steps.back());


      ABCRL<G, TotalRewardStatistic<Vector, int>, M, AbstractPolicy<Vector, int>, Vector, int> abcrl;
      std::vector<M> samples;
      std::vector<real> values;
      abcrl.GetSample(training_data_list,
		      policy_list,
		      options.gamma,
		      options.epsilon,
		      options.delta,
		      options.R_max,
		      options.n_trajectories,
		      options.n_samples,
		      samples,
		      values);
      Options estimation_options = options;
      estimation_options.n_training = options.n_model_rollouts;
      AbstractPolicy<Vector, int>* sampled_policy
	=  getLSPIPolicy(&samples[0],
			 *policy_list[0],
			 NULL,
			 RBFs,
			 estimation_options);
      policy_list.push_back(sampled_policy);
    } else {
      training_rollouts.StartingDistributionSampling(n_episodes_per_iteration, rollout_horizon);
      printf("%f %f %d # total, discounted reward, steps LSPI\n",
	     training_rollouts.total_rewards.back(),
	     training_rollouts.discounted_rewards.back(),
	     training_rollouts.total_steps.back());

      AbstractPolicy<Vector, int>* lspi_policy
	=  getLSPIPolicy(NULL,
			 policy,
			 &training_rollouts,
			 RBFs,
			 options);
      epsilon_greedy *= pow(0.99, (real) training_rollouts.total_steps.back());
      lspi_policy->setEpsilonGreedy(epsilon_greedy);
      logmsg("New epsilon %f\n", epsilon_greedy);
      policy_list.push_back(lspi_policy);
      training_rollouts.policy = lspi_policy;
    }
  }
  for (uint i=0; i<policy_list.size(); ++i) {
    delete policy_list[i];
  }
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
    --environment:           {MountainCar, Pendulum, Puddle, Bicycle, CartPole, Acrobot}\n\
    --discount:              reward discounting in [0,1]\n\
    --threshold:             statistic threshold\n\
    --n_trajectories:        number of trajectories per sample\n\
    --n_samples:             number of sampled models\n\
    --n_training:            number of training trajectories\n\
    --n_testing:             number of test trajectories\n\
    --seed:                  seed all the RNGs with this\n\
    --seed_file:             select a binary file to choose seeds from (use in conjunction with --seed to select the n-th seed in the file)\n\
    --grid:                  number of grid intervals for LSTD\n\
    --scale:                 RBF scale for LSTD\n\
    --n_evaluations:         number of evaluations to run\n\
    --reuse_training_data:   reuse the training data in LSPI samples\n\
    --online:                do the online test\n\
    --sampling:              use sampling\n\
    --delta:                 error bound probability\n\
    --Rmax:                  maximum reward value\n\
	--c_uct;                 uct factor\n\
	--n_rollouts_uct;        number of rollouts to perform with UCT\n\
	--lambda_uct;            lambda mixing for uct\n\
	--stepsize_uct;          step size\n\
	--depth_uct;             maximum tree depth\n\
	--grid_uct;              grid utc\n\
\n";

int main(int argc, char* argv[])
{
  ulong seed = time(NULL);
  char* seed_filename = 0;

  bool online_test = false;


  MersenneTwisterRNG rng;
  Options options(rng);
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
	{"sampling", no_argument, 0, 0}, //15
	{"seed_file", required_argument, 0, 0}, //16
	{"delta", required_argument, 0, 0}, //17
	{"Rmax", required_argument, 0, 0}, //18
	{"c_uct", required_argument, 0, 0}, //19
	{"n_rollouts_uct", required_argument, 0, 0}, //20
	{"lambda_uct", required_argument, 0, 0}, //21
	{"stepsize_uct", required_argument, 0, 0}, //22
	{"depth_uct", required_argument, 0, 0}, //23
	{"grid_uct", required_argument, 0, 0}, //24
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
	case 15: options.sampling = true; break;
	case 16: seed_filename = optarg; break;
	case 17: options.delta = atof(optarg); break;
	case 18: options.R_max = atof(optarg); break;
	case 19: options.c_uct = atof(optarg); break;
	case 20: options.n_rollouts_uct = atoi(optarg); break;
	case 21: options.lambda_uct = atof(optarg); break;
	case 22: options.stepsize_uct = atof(optarg); break;
	case 23: options.depth_uct = atoi(optarg); break;
	case 24: options.grid_uct = atoi(optarg); break;
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

  options.ShowOptions();

  if (online_test) {
    for (int i=0; i<options.n_evaluations; ++i) {
      logmsg("run %d/%d\n", i, options.n_evaluations);
      if (!strcmp(options.environment_name, "MountainCar")) {
	logmsg("Testing mountain car\n");
	RunOnlineTest<MountainCarGenerator, MountainCar>(options);
      } else if (!strcmp(options.environment_name, "Pendulum")) {
	logmsg("Testing pendulum\n");
	RunOnlineTest<PendulumGenerator, Pendulum>(options);
      } else if (!strcmp(options.environment_name, "PuddleWorld")) {
	logmsg("Testing puddle world\n");
	RunOnlineTest<PuddleWorldGenerator, PuddleWorld>(options);
      } else if (!strcmp(options.environment_name, "Bicycle")) {
	logmsg("Testing Bike\n");
	RunOnlineTest<BikeGenerator, Bike>(options);
      } else if (!strcmp(options.environment_name, "CartPole")) {
	logmsg("Testing CartPole\n");
	RunOnlineTest<CartPoleGenerator, CartPole>(options);
      } else if (!strcmp(options.environment_name, "Acrobot")) {
	logmsg("Testing Acrobot\n");
	RunOnlineTest<AcrobotGenerator, Acrobot>(options);
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
      } else if (!strcmp(options.environment_name, "Bicycle")) {
	logmsg("Testing Bike\n");
	RunTest<BikeGenerator, Bike>(options);
      } else if (!strcmp(options.environment_name, "CartPole")) {
	logmsg("Testing CartPole\n");
	RunTest<CartPoleGenerator, CartPole>(options);
      } else if (!strcmp(options.environment_name, "Acrobot")) {
	logmsg("Testing Acrobot\n");
	RunTest<AcrobotGenerator, Acrobot>(options);
      } else {
	fprintf(stderr, "Invalid environment name %s\n", options.environment_name);
	exit(-1);
      }
    }
  }
  return 0;
}
#endif
