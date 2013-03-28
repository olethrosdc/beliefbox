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

#ifndef GRADIENT_BRL_H
#define GRADIENT_BRL_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "MDPModel.h"
#include "ValueIteration.h"
#include <vector>

/// \ingroup ReinforcementLearning
/// @{
    
/** Gradient Bayesian Reinforcement Learning.

    The idea is that we sample an MDP at a time, and then move our
    value function closer to the value function of that MDP.  Thus,
    our value function moves to the expected value function!

    Unfortunately, calculating the value of an MDP is slow, so just
    perform a single TD update.

 */
class GradientBRL : public OnlineAlgorithm<int, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real epsilon; ///< randomness
    real step_size; ///< step size for the gradient
    int current_state; ///< current state
    int current_action; ///< current action
    MDPModel* model; ///< pointer to the base MDP model
    ValueIteration* value_iteration; ///< only one value iteration required. 
    std::vector<real> tmpQ;
    RandomNumberGenerator* rng; ///< random number generator to draw samples from
    int T; ///< time passed
    int update_interval; ///< update interval for policy
    int next_update; ///< next time at which we use a new policy
    bool use_upper_bound; ///< use upper bounds to take actions if true
    bool use_sampling_threshold; ///< use a threshold for resampling
    real sampling_threshold; ///< value of the threshold

public:
    const DiscreteMDP* mdp_sample; ///< sampled model
    Vector weights; ///< probability vector of MDPs
    GradientBRL(int n_states_,
                int n_actions_,
                real gamma_,
                real epsilon_,
                real step_size_,
                MDPModel* model_,
                RandomNumberGenerator* rng_,
                bool use_upper_bound_ = false);
    virtual ~GradientBRL();
    virtual void Reset();
    /// Full observation
    virtual real Observe (int state, int action, real reward, int next_state, int next_action);
    /// Partial observation 
    virtual real Observe (real reward, int next_state, int next_action);

    /// Sample a new set of MDPs
    void Resample();

    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, int next_state);


    virtual real getValue (int state, int action)
    {
        return value_iteration->getValue(state, action);
    }
    
    /** Set the rewards to Singular distributions.

        Since this is a Bayesian approach, we can simply set the belief about the reward in each state to be a singular distribution.
    */
    virtual void setFixedRewards(const Matrix& rewards)
    {
        model->setFixedRewards(rewards);
    }

};


/// @}
#endif

