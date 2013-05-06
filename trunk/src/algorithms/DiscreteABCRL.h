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

#ifndef DISCRETE_ABCRL_H
#define DISCRETE_ABCRL_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "Demonstrations.h"
#include "MultiMDPValueIteration.h"
#include "ValueIteration.h"
#include "Environment.h"
#include <vector>

/// \ingroup ReinforcementLearning
/// @{
    
/** Direct model-based reinforcement learning.
    
    This class maintains a set of (discrete) MDPs.
    It selects actions based on a set of _sampled_ MDPs.
    
    The sampling is performed using ABC, using a variant of the
    algorithm described in "ABC Reinforcement Learning" Christos
    Dimitrakakis and Nikolaos Tziortziotis ICML 2013.

    The lower bound algorithm is described in the paper:
    "Robust Bayesian Reinforcement Learning through Tight Lower Bounds",
    Christos Dimitrakakis
    EWRL 2011.
*/
class DiscreteABCRL : public OnlineAlgorithm<int, int>
{
protected:
  const int n_states; ///< number of states
  const int n_actions; ///< number 
  real gamma; ///< discount factor
  real epsilon; ///< randomness
  int current_state; ///< current state
  int current_action; ///< current action
  EnvironmentGenerator<int, int>* generator; ///< generator
  Demonstrations<int, int> demonstrations; ///< demonstrations
  std::vector<DiscretePolicy*> policies;
  std::vector<ValueIteration*> value_iteration; ///< value iteration on each separate model
  MultiMDPValueIteration* multi_value_iteration; ///< multi-MDP value iteration
  std::vector<real> tmpQ;
  Vector VU; ///< upper bound value
  Vector VL; ///< lower bound value
  Matrix QU; ///< upper bound Q-value
  Matrix QL; ///< lower bound Q-value
  int max_samples; ///< maximum number of samples to take
  RandomNumberGenerator* rng; ///< random number generator to draw samples from
  int T; ///< time passed
  int update_interval; ///< update interval for policy
  int next_update; ///< next time at which we use a new policy
  bool use_upper_bound; ///< use upper bounds to take actions if true
  bool use_sampling_threshold; ///< use a threshold for resampling
  real sampling_threshold; ///< value of the threshold

public:
  std::vector<const DiscreteMDP*> mdp_list; ///< list of sampled models
  Vector weights; ///< probability vector of MDPs
  DiscreteABCRL(int n_states_,
                int n_actions_,
                real gamma_,
                real epsilon_,
                EnvironmentGenerator<int, int>* generator_,
                RandomNumberGenerator* rng_,
                int max_samples_ = 1,
                bool use_upper_bound_ = false);
  virtual ~DiscreteABCRL();
  virtual void Reset();
  /// Full observation
  virtual real Observe (int state, int action, real reward, int next_state, int next_action);
  /// Partial observation 
  virtual real Observe (real reward, int next_state, int next_action);
  DiscreteMDP* GenerateMDP() const;
  /// Sample a new set of MDPs
  void Resample();

  /// Get an action using the current exploration policy.
  /// it calls Observe as a side-effect.
  virtual int Act(real reward, int next_state);


  virtual real getValue (int state, int action)
  {
    if (use_upper_bound) {
      return UpperBound(state, action);
    } else {
      return LowerBound(state, action);
    }
  }
    
  void CalculateUpperBound(real accuracy, int iterations);

  void CalculateLowerBound(real accuracy, int iterations);

  inline real UpperBound(int state)
  {
    Vector Q(n_actions);
    for (int i=0; i<n_actions; ++i) {
      Q(i) = UpperBound(state, i);
    }
    return Max(Q);
  }

  inline real LowerBound(int state)
  {
    Vector Q(n_actions);
    for (int i=0; i<n_actions; ++i) {
      Q(i) = LowerBound(state, i);
    }
    return Max(Q);
  }

  inline real UpperBound(int state, int action)
  {
    return QU(state, action);
  }

  inline real LowerBound(int state, int action)
  {
    return QL(state, action);
  }
    
  /** Set the rewards to Singular distributions.

      Since this is a Bayesian approach, we can simply set the belief about the reward in each state to be a singular distribution.
  */
  virtual void setFixedRewards(const Matrix& rewards)
  {
    //model->setFixedRewards(rewards);
#if 0
    logmsg("Setting reward matrix\n");
    rewards.print(stdout);
    model->ShowModel();
#endif
  }

  virtual void setSamplingThreshold(real sampling_threshold_)
  {
    use_sampling_threshold = true;
    sampling_threshold = sampling_threshold_;
    assert(sampling_threshold >= 0.0 && sampling_threshold <= 1.0);
  }
};


/// @}
#endif

