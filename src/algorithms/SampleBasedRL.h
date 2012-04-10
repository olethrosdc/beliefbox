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
    
/** Direct model-based reinforcement learning.
    
    This class maintains a model of the (discrete) MDP.
    It selects actions based on a set of _sampled_ MDPs.
    
    This is described in the paper

    "Robust Bayesian Reinforcement Learning through Tight Lower Boudns"
    EWRL 2012.
*/
class SampleBasedRL : public OnlineAlgorithm<int, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real epsilon; ///< randomness
    int state; ///< current state
    int action; ///< current action
    MDPModel* model;
    std::vector<const DiscreteMDP*> mdp_list;
    std::vector<ValueIteration*> value_iteration;
    MultiMDPValueIteration* multi_value_iteration;
    std::vector<real> tmpQ;
    int max_samples;
    RandomNumberGenerator* rng;
    int T;
    int update_interval;
    int next_update;
    bool use_upper_bound;
public:
    SampleBasedRL(int n_states_,
                  int n_actions_,
                  real gamma_,
                  real epsilon_,
                  MDPModel* model_,
                  RandomNumberGenerator* rng_,
                  int max_samples_ = 1,
                  bool use_upper_bound_ = false);
    virtual ~SampleBasedRL();
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
        if (use_upper_bound) {
            return UpperBound(state, action);
        } else {
            return LowerBound(state, action);
        }
    }

    inline real UpperBound(int state, int action)
    {
        real Q = 0.0;
        for (int i=0; i<max_samples; i++) {
            Q += value_iteration[i]->getValue(state, action);
        }
        return Q / (real) max_samples;
    }

    inline real LowerBound(int state, int action)
    {
        return multi_value_iteration->getValue(state, action);
    }
    
};


/// @}
#endif

