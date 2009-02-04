// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef Q_LEARNING_H
#define Q_LEARNING_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "real.h"
#include "ExplorationPolicy.h"
#include "OnlineAlgorithm.h"
#include <vector>

class QLearning : public OnlineAlgorithm<int,int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real lambda; ///< eligibility trace decay rate
    real alpha; ///< learning rate 
    VFExplorationPolicy* exploration_policy; ///< exploration policy
    real initial_value; ///< initial value for Q values
    real baseline; ///< baseline reward

    Matrix Q;
    Matrix el;

    int state; ///< current state
    int action; ///< current action
public:
    QLearning(int n_states_,
              int n_actions_,
              real gamma_,
              real lambda_,
              real alpha_,
              VFExplorationPolicy* exploration_policy_,
              real initial_value_= 0.0,
              real baseline_ = 0.0);
    virtual ~QLearning()
    {
    }
    virtual void Reset();
    virtual real Observe (int action, int next_state, real reward);
    virtual real Observe (real reward, int next_state, int next_action);
    virtual int Act(real reward, int next_state);
    virtual real getValue (int s, int a)
    {
        return Q(state, action);
    }

    
};

#endif

