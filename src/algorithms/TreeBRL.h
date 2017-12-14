// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SAMPLE_BASED_RL_H
#define SAMPLE_BASED_RL_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "MDPModel.h"
#include "MultiMDPValueIteration.h"
#include "ValueIteration.h"
#include <vector>

/// \ingroup ReinforcementLearning
/// @{
    
/** Direct model-based reinforcement learning using trees.
  
*/
class TreeBRL : public OnlineAlgorithm<int, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real epsilon; ///< randomness
    int current_state; ///< current state
    int current_action; ///< current action
    MDPModel* model; ///< pointer to the base MDP model
    int horizon; ///< maximum number of samples to take
    RandomNumberGenerator* rng; ///< random number generator to draw samples from
    int T; ///< time passed
public:
  class BeliefState
  {
  protected:
    TreeBRL& tree;
    Matrix P;
    int state;
  public:
    BeliefState(TreeBRL& tree, MDPModel* model, int state);
    void ExpandAllActions()
    {
      const DiscreteMDP* mdp = model->generate();
      for (int a=0; a<tree.n_actions; ++a) {
	for (int s=0; s<tree.n_states; ++s) {
	  x
	}
      }
    }
    
  };
  TreeBRL(int n_states_,
	  int n_actions_,
	  real gamma_,
	  MDPModel* model_,
	  RandomNumberGenerator* rng_,
	  int horizon_ = 1,
	  bool use_upper_bound_ = false);
    virtual ~TreeBRL();
    virtual void Reset();
    /// Full observation
    virtual real Observe (int state, int action, real reward, int next_state, int next_action);
    /// Partial observation 
    virtual real Observe (real reward, int next_state, int next_action);
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, int next_state);
    virtual real getValue (int state, int action)
    {
      return BQ(state, action);
    }
    /** Set the rewards to Singular distributions.

        Since this is a Bayesian approach, we can simply set the belief about the reward in each state to be a singular distribution.
     */
    virtual void setFixedRewards(const Matrix& rewards)
    {
        model->setFixedRewards(rewards);
#if 0
        logmsg("Setting reward matrix\n");
        rewards.print(stdout);
		model->ShowModel();
#endif
    }
};


/// @}
#endif

