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

#ifndef TD_MB_H
#define	TD_MB_H


#include <map>
#include <string>
#include <utility>
#include <vector>

#include "limits.h"
#include "real.h"

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "NormalDistribution.h"
#include "OnlineAlgorithm.h"


using namespace std;


class TdMb : public OnlineAlgorithm<int, int>
{

protected:
	const int n_states; ///< number of states
	const int n_actions; ///< number 

	real alpha; ///< learning rate 
	real gamma; ///< discount factor
	real lambda; ///< eligibility trace decay rate

	VFExplorationPolicy* exploration_policy; ///< exploration policy
	real initial_value; ///< initial value for Q values
	real baseline; ///< baseline reward

	Matrix Q; ///< values of all states and actions

	Matrix E; ///< eligibility traces for all states and actions
	Matrix LAMBDA; ///< variable lambda
	map < pair<int, int>, vector<int> > TRIALS; ///< counts the number of outcomes
	Matrix N; ///< counts the number of trials
	real lambdaParameter;

	int t, T;
	vector< pair< pair< int, int >, double > > episode; ///< vector of ((state, action), reward)

	int state; ///< current state
	int action; ///< current action


public:
	TdMb(int n_states_,
		int n_actions_,
		real gamma_,
		real lambda_,
		real alpha_,
		VFExplorationPolicy* exploration_policy_,
		real initial_value_ = 0.0,
		real baseline_ = 0.0);
	virtual ~TdMb();

	virtual void Reset();

	/// Full TD-MB observation (no eligibility traces)
	virtual real Observe(int state, int action, real reward, int next_state, int next_action);

	/// Partial TD-MB observation (can be used with eligibility traces)
	virtual real Observe(real reward, int next_state, int next_action);

	/// Get an action using the current exploration policy.
	/// It calls Observe as a side-effect.
	virtual int Act(real reward, int next_state);

	virtual real getValue(int state, int action)
	{
		return Q(state, action);
	}
};


#endif	/* TD_MB_H */
