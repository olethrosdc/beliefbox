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
#include "NormalDistribution.h"
#include "OnlineAlgorithm.h"


using namespace std;


class TdBma : public OnlineAlgorithm<int, int>
{

protected:
	const int N_STATES; ///< number of states
	const int N_ACTIONS; ///< number 

	const real ALPHA; ///< learning rate 
	const real GAMMA; ///< discount factor
	real lambda; ///< eligibility trace decay rate

	VFExplorationPolicy* exploration_policy; ///< exploration policy
	const real INITIAL_VALUE; ///< initial value for Q values
	const real BASELINE; ///< baseline reward

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
	TdBma(const int n_states_,
		const int n_actions_,
		const real gamma_,
		const real lambda_,
		const real alpha_,
		VFExplorationPolicy* exploration_policy_,
		const real initial_value_ = 0.0,
		const real baseline_ = 0.0);
	virtual ~TdBma();

	virtual void Reset();

	/// Full TD-BMA observation (no eligibility traces)
	virtual real Observe(const int state, const int action, const real reward, const int next_state, const int next_action);

	/// Partial TD-BMA observation (can be used with eligibility traces)
	virtual real Observe(const real reward, const int next_state, const int next_action);

	/// Get an action using the current exploration policy.
	/// It calls Observe as a side-effect.
	virtual int Act(const real reward, const int next_state);

	virtual real getValue(const int state, const int action)
	{
		return Q(state, action);
	}
};


#endif	/* TD_BMA_H */
