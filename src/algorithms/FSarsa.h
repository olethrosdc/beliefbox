// -*- Mode: c++ -*-
// copyright (c) 2021 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FSARSA_H
#define FSARSA_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include <vector>

/// Sarsa with features
class Sarsa : public OnlineAlgorithm<int, Vector>
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

	Vector V;
    Matrix Q;
    Matrix el;

    int state; ///< current state
    int action; ///< current action

public:
    Sarsa(int n_states_,
          int n_actions_,
          real gamma_,
		  BasisSet& basis,
          real lambda_=0.0,
          real alpha_=0.5,
          VFExplorationPolicy* exploration_policy_=NULL,
          real initial_value_= 0.0,
          real baseline_ = 0.0);
    virtual ~Sarsa();
    virtual void Reset();
    /// Full SARSA observation (no eligibility traces)
    virtual real Observe (const Vector& state, const int& action, real reward, const Vector& next_state, const int& next_action);
    /// Partial SARSA observation (can be used with eligibility traces)
    virtual real Observe (real reward, const Vector& next_state, const int& next_action);
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, const Vector& next_state);

    virtual real getValue (const Vector& state, const int& action)
    {
        return Q(state, action);
    }
	virtual real getValue (int state)
    {
        return V(state);
    }

    
};

#endif

