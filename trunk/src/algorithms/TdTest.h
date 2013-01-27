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

#ifndef TD_TEST_H
#define	TD_TEST_H


#include <map>
#include <string>
#include <utility>
#include <vector>

#include "limits.h"
#include "real.h"
#include "SmartAssert.h"

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "NormalDistribution.h"
#include "OnlineAlgorithm.h"


using namespace std;


class TdTest : public OnlineAlgorithm<int, int>
{

protected:
	const int N_STATES; ///< number of states
	const int N_ACTIONS; ///< number of actions

	const real ALPHA; ///< learning rate 
	const real GAMMA; ///< discount factor
	const real DEFAULT_LAMBDA; ///< default eligibility trace decay rate

	VFExplorationPolicy* exploration_policy; ///< exploration policy
	const real INITIAL_VALUE; ///< initial value for Q values
	const real BASELINE; ///< baseline reward

	Matrix Q; ///< contains the values of every state-action pair

	Matrix E; ///< eligibility traces for all states and actions
	Matrix LAMBDA; ///< eligibility trace decay rate for a state and an action
	Matrix N; ///< counts the number of trials
	const real U; ///< Dirichlet parameter
	map < pair<int, int>, vector<int> > TRIALS; ///< counts the number of outcomes
	Matrix MEANS; ///< stores the mean value of every state-action pair

	int state; ///< current state
	int action; ///< current action

private:

	/// Get the probability from the Dirichlet distribution

	double GetProbability(int i, int n)
	{
		assert(n != 0 && n >= i);
		return ((double) (i + U)) / ((double) (n * (1 + U)));
	}

	/// Get the maximum value of Q(s, a) for state s

	double QMax(int s)
	{
		double qMax = 0;
		for(int a = 0; a < N_ACTIONS; ++a)
		{
			if(qMax < Q(s, a))
			{
				qMax = Q(s, a);
			}
		}
		return qMax;
	}

public:
	TdTest(const int n_states_,
		const int n_actions_,
		const real gamma_,
		const real lambda_,
		const real alpha_,
		VFExplorationPolicy * const exploration_policy_,
		const real initial_value_ = 0.0,
		const real baseline_ = 0.0);
	virtual ~TdTest();

	virtual void Reset();

	/// Full TD-MB observation
	virtual real Observe(const int state, const int action, const real reward, const int next_state, const int next_action);

	/// Partial TD-MB observation (can be used with eligibility traces)
	virtual real Observe(const real reward, const int next_state, const int next_action);

	/// Get an action using the current exploration policy.
	/// It calls Observe as a side-effect.
	virtual int Act(const real reward, const int next_state);

	virtual real getValue(const int state, const int action)
	{
		return Q(state, action);
	}
};


#endif	/* TD_TEST_H */
