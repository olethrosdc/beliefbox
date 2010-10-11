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
#include <vector>

/// \ingroup ReinforcementLearning
/// @{
	
/** Direct model-based reinforcement learning.
	
	This class maintains a model of the (discrete) MDP.
    It selects actions based on a set of _sampled_ MDPs.
	
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
    MultiMDPValueIteration* value_iteration;
    std::vector<real> tmpQ;
    int max_samples;
    int T;
    int update_interval;
public:
    SampleBasedRL(int n_states_,
                 int n_actions_,
                 real gamma_,
                 real epsilon_,
                 MDPModel* model_,
                 int max_samples_ = 1);
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
        return value_iteration->getValue(state, action);
    }
    
};


/// @}
#endif

