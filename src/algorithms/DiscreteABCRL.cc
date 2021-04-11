// -*- Mode: c++ -*-
// copyright (c) 2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteABCRL.h"

DiscreteABCRL::DiscreteABCRL(int n_states_,
                             int n_actions_,
                             real gamma_,
                             real epsilon_,
                             EnvironmentGenerator<int, int>* generator_,
                             RandomNumberGenerator* rng_,
                             int max_samples_,
                             int n_iterations_,
                             bool use_upper_bound_)
  : n_states(n_states_),
    n_actions(n_actions_),
    gamma(gamma_),
    epsilon(epsilon_),
    generator(generator_),
    demonstrations(false),
    tmpQ(n_actions),
    VU(n_states),
    VL(n_states),
    QU(n_states, n_actions),
    QL(n_states, n_actions),
    max_samples(max_samples_),
    rng(rng_),
    T(0),
    update_interval(n_states + 1),
    next_update(0),
    use_upper_bound(use_upper_bound_),
    use_sampling_threshold(false),
    sampling_threshold(0.1),
    n_iterations(n_iterations_),
    weights(max_samples)
{
  printf("# Starting Discrete-ABC-RL with %d samples, update interval %d\n",
         max_samples, update_interval);
  
  current_state = -1;
  
  real w_i = 1.0 / (real) max_samples;
  mdp_list.resize(max_samples);
  value_iteration.resize(max_samples);
  printf("# Generating mean MDP\n");
  //mdp_list[0] = model->getMeanMDP();
  for (int i=0; i<max_samples; ++i) {
    printf("# Generating sampled MDP\n");
    mdp_list[i] = GenerateMDP();
    weights[i] = w_i;
    value_iteration[i] = new ValueIteration(mdp_list[i], gamma);
  }

  printf ("# Setting up MultiMPDValueIteration\n");
  multi_value_iteration = new MultiMDPValueIteration(weights, mdp_list, gamma);
  printf ("# Testing MultiMPDValueIteration\n");
  multi_value_iteration->ComputeStateActionValues(0,1);
}
DiscreteABCRL::~DiscreteABCRL()
{

  for (int i=0; i<max_samples; ++i) {
    delete mdp_list[i];
    delete value_iteration[i];
  }
  delete multi_value_iteration;
}

void DiscreteABCRL::Reset()
{
  current_state = -1;
  next_update = T;
  Resample();
  //model->Reset();
}
/// Full observation
real DiscreteABCRL::Observe (const int& state, const int& action, real reward, const int& next_state, const int& next_action)
{
  if (state >= 0) {
    demonstrations.Observe(state, action, reward);
  } else {
    demonstrations.NewEpisode();
  }
  current_state = next_state;
  current_action = next_action;
  return 0.0;
}

/// Partial observation 
real DiscreteABCRL::Observe (real reward, const int& next_state, const int& next_action)
{
  if (current_state < 0) {
    demonstrations.NewEpisode();
  } else {
    demonstrations.Observe(current_state, current_action, reward);
  }
  current_state = next_state;
  current_action = next_action;
  return 0.0;
}

DiscreteMDP* DiscreteABCRL::GenerateMDP() const
 {
  int iter = 0;
  real min_error = INF;
  DiscreteMDP* mdp = NULL;
  while (iter < n_iterations) {
    DiscreteEnvironment* environment = generator->Generate();
    Demonstrations<int, int> test_demos(false);
    if (!mdp) {
      mdp = environment->getMDP();
    }
    for (uint i=0; i<policies.size(); ++i) {
      test_demos.Simulate(*environment, *policies[i], gamma, demonstrations.length(i));
    }
    real mean = 0;
    for (uint i=0; i<policies.size(); ++i)  {
      printf("D %f\n", demonstrations.discounted_reward(i));
      mean += demonstrations.discounted_reward(i);
    }
    mean /= (real) demonstrations.size();
    real test_mean = 0;
    for (uint i=0; i<test_demos.size(); ++i)  {
      printf("T %f\n", test_demos.discounted_reward(i));
      test_mean += test_demos.discounted_reward(i);
    }
    test_mean /= (real) test_demos.size();
    ++iter;
    real error = fabs(mean - test_mean);
    logmsg("err: %f (%f %f)\n", error, mean, test_mean);
    if (error < min_error) {
      min_error = error;
      delete mdp;
      mdp = environment->getMDP();
    }
    delete environment;
  }
  
  logmsg("utility error: %f\n", min_error);
  return mdp;

}

void DiscreteABCRL::Resample()
{
  for (int i=0; i<max_samples; ++i) {
    delete mdp_list[i];
    mdp_list[i] = GenerateMDP();
  }
}


void DiscreteABCRL::CalculateUpperBound(real accuracy, int iterations)
{
  for (int j=0; j<max_samples; ++j) {
    value_iteration[j]->setMDP(mdp_list[j]);
    value_iteration[j]->ComputeStateValuesStandard(accuracy, iterations);
  }
    
  real Z = 1.0 / (real) max_samples;
  for (int s=0; s<n_states; ++s) {
    for (int a=0; a<n_actions; ++a) {
      QU(s,a) = 0;
      for (int i=0; i<max_samples; i++) {
        QU(s, a) += value_iteration[i]->getValue(s, a);
      }
    }
  }
  QU *= Z;
  for (int s=0; s<n_states; ++s) {
    VU(s) = QU(s, 0);
    for (int a=1; a<n_actions; ++a) {
      VU(s) = std::max<real>(VU(s), QU(s, a));
    }
  }
}

void DiscreteABCRL::CalculateLowerBound(real accuracy, int iterations)
{
  multi_value_iteration->setMDPList(mdp_list);
  multi_value_iteration->ComputeStateValues(accuracy, iterations);
    
  for (int s=0; s<n_states; ++s) {
    for (int a=0; a<n_actions; ++a) {
      QL(s,a) = multi_value_iteration->getValue(s, a);
    }
  }
  for (int s=0; s<n_states; ++s) {
    VL(s) = QL(s, 0);
    for (int a=1; a<n_actions; ++a) {
      VL(s) = std::max<real>(VL(s), QL(s, a));
    }
  }
}


/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int DiscreteABCRL::Act(real reward, const int& next_state)
{
  assert(next_state >= 0 && next_state < n_states);
  T++;

  // update the model
  if (current_state >= 0) {
    demonstrations.Observe(current_state, current_action, reward);
  } else {
    demonstrations.NewEpisode();
  }
  current_state = next_state;

  // Update MDPs
  //mdp_list[0] = model->getMeanMDP();
  // Do note waste much time generating MDPs
    
  bool do_update = false;
  if (use_sampling_threshold) {
    for (int i=0; i<max_samples; ++i) {
      real p = mdp_list[i]->getTransitionProbability(current_state, current_action, next_state);
      weights[i] *= p;
    }
    weights /= weights.Sum();
    if (Max(weights) > 1.0 - sampling_threshold) {
      //printf("Minimum: %f\n", Min(weights));
      //weights.print(stdout);
      do_update = true;
      for (int i=0; i<max_samples; ++i) {
        weights[i] = 1.0 / (real) max_samples;
      }
    }
  } else {
    if (T >= next_update) {    
      do_update = true;
    }
  }
  if (do_update) {    
    //printf("# update: %d\n", T);
    //model->ShowModel();
    update_interval += 1;//(int) (ceil)(1.01*(double) T);
    next_update = T + update_interval;
    Resample();
    if (use_upper_bound) {
      CalculateUpperBound(1e-3, 1e3);
    } else {
      CalculateLowerBound(1e-3, 1e3);
    }
  }

  // update values    
  if (use_upper_bound) {
    //CalculateUpperBound(0, 1);
    for (int i=0; i<n_actions; i++) {
      tmpQ[i] = UpperBound(next_state, i);
    }
  } else {
    //CalculateLowerBound(0, 1);
    for (int i=0; i<n_actions; i++) {
      tmpQ[i] = LowerBound(next_state, i);
    }

  }

  int next_action;
  real epsilon_t = epsilon / (1.0 + sqrt((real) T));

  // choose action
  if (urandom()<epsilon_t) {
    next_action = rng->discrete_uniform(n_actions);
    //printf("RANDOM %d\n", next_action);
  } else {
    next_action = ArgMax(tmpQ);
  }
  current_action = next_action;
  //printf("%f %d #epsilon\n", epsilon_t, action);
    
  return current_action;
}


