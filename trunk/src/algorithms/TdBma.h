/* -*- Mode: C++; -*- */
/* (c) 2012 Florian Barras */
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TDA_BMA_H
#define TDA_BMA_H


#include <string>
#include <utility>
#include <vector>

#include "real.h"

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "OnlineAlgorithm.h"


using namespace std;

class TdBma : public OnlineAlgorithm<int, int>
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

	Matrix N;
	Matrix Means;
	Matrix Variances;
	Matrix Deviations;

	int t, T;
	vector< pair< pair< int, int >, double > > episode; ///< vector of ((state, action), reward)

	int state; ///< current state
	int action; ///< current action


public:
	TdBma(int n_states_,
			int n_actions_,
			real gamma_,
			real lambda_,
			real alpha_,
			VFExplorationPolicy* exploration_policy_,
			real initial_value_ = 0.0,
			real baseline_ = 0.0);
	virtual ~TdBma();

	virtual void Reset();

	/// Full TD-BMA observation (no eligibility traces)
	virtual real Observe(int state, int action, real reward, int next_state, int next_action);

	/// Partial TD-BMA observation (can be used with eligibility traces)
	virtual real Observe(real reward, int next_state, int next_action);

	/// Get an action using the current exploration policy.
	/// it calls Observe as a side-effect.
	virtual int Act(real reward, int next_state);

	virtual real getValue(int state, int action) {
		return Q(state, action);
	}
};


#endif
