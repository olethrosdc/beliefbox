// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef UCRL2_H
#define UCRL2_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "MDPModel.h"
#include "OptimisticValueIteration.h"
#include "RandomNumberGenerator.h"
#include <vector>

/** 
    \ingroup ReinforcementLearning
*/
/*@{*/
	
/** Direct model-based reinforcement learning.
	
	This class maintains a model of the (discrete) MDP.
	
*/
class UCRL2 : public OnlineAlgorithm<int, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real confidence_interval; ///< confidence interval
    int state; ///< current state
    int action; ///< current action
    DiscreteMDPCounts* model; ///< stores what is known about the model
    OptimisticValueIteration* value_iteration;
    std::vector<real> tmpQ;
	RandomNumberGenerator* rng;
    int total_steps;
    int update_interval;
    int next_update;
    Matrix rewards;
    bool known_rewards;
    int n_resets;
public:
    UCRL2(int n_states_,
          int n_actions_,
          real gamma_,
          DiscreteMDPCounts* model_,
          RandomNumberGenerator* rng_,
		  real delta);
    virtual ~UCRL2();
    virtual void Reset();
    /// Full observation
    virtual real Observe (const int& state, const int& action, real reward, const int& next_state, const int& next_action);
    /// Partial observation 
    virtual real Observe (real reward, const int& next_state, const int& next_action);
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, const int& next_state);

    virtual real getValue (const int& state, const int& action)
    {
        return value_iteration->getValue(state, action);
    }

    virtual void setFixedRewards(const Matrix& rewards)
    {
        known_rewards = true;
        this->rewards = rewards;
		model->setFixedRewards(rewards);
#if 0
		rewards.print(stdout);
		model->ShowModel();
#endif
		value_iteration->ComputeStateValuesKnownRewards(confidence_interval, 1e-6, -1);
    }
    
};


/*@}*/
#endif

