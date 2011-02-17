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

/** A simple implementation of \f$Q(\lambda)\f$.
	
	This is an online algorithm operating on discrete state-action
	spaces. It is an off-policy algorithm.

	The main additional parameter is the exploration policy, which
	can be defined separately from the learning algorithm itself.

	@see Sutton and Barto 1998: "Introduction to reinforcement learning".
 */
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

    Matrix Q; ///< The matrix of Q values
    Matrix el; ///< The matrix of eligibility traces

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
	/// Destructor
    virtual ~QLearning()
    {
    }
    virtual void Reset();
    virtual real Observe (int action, int next_state, real reward);
    virtual real Observe (real reward, int next_state, int next_action);
    virtual int Act(real reward, int next_state);
	/// Get value of state-action
	virtual real getValue (int s, int a)
    {
        return Q(s, a);
    }
	virtual real& QValue (int s, int a)
    {
        return Q(s, a);
    }
    Matrix getQMatrix()
    {
        return Q;
    }
};

#endif

