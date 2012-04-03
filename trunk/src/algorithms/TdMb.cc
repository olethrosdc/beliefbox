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

#include "TdMb.h"

TdMb::TdMb(int n_states_,
	int n_actions_,
	real gamma_,
	real lambda_,
	real alpha_,
	VFExplorationPolicy* exploration_policy_,
	real initial_value_,
	real baseline_)
	: n_states(n_states_),
	n_actions(n_actions_),
	gamma(gamma_),
	lambda(lambda_),
	alpha(alpha_),
	exploration_policy(exploration_policy_),
	initial_value(initial_value_),
	baseline(baseline_),
	Q(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	E(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	LAMBDA(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	N(n_states_, n_actions_, Matrix::CHECK_BOUNDS)
{
	assert(alpha >= 0 && alpha <= 1);
	assert(gamma >= 0 && gamma <= 1);
	assert(lambda >= 0 && lambda <= 1);

	lambdaParameter = 0.5;

	Reset();

	exploration_policy->setValueMatrix(&Q);
}

TdMb::~TdMb() { }

void TdMb::Reset()
{
	state = -1;
	action = -1;

	T = n_states;
	t = 0;

	//episode();

	for(int s = 0; s < n_states; ++s)
	{
		for(int a = 0; a < n_actions; ++a)
		{
			Q(s, a) = initial_value;

			E(s, a) = 0.0;
			LAMBDA(s, a) = lambda;

			vector<int> sa_trials(n_states, 0);
			pair<int, int> sa(s, a);
			TRIALS[sa] = sa_trials;

			N(s, a) = 0;
		}
	}

}

/// Full TD-MB observation (no eligibility traces)

real TdMb::Observe(int state, int action, real reward, int next_state, int next_action)
{
	// o Compututation

	// - Delta
	real delta = reward + gamma * Q(next_state, next_action) - Q(state, action);

	// - Lambda
	pair<int, int> sa(state, action);
	vector<int> trials = TRIALS[sa];

	trials[next_state] += 1;
	N(state, action) += 1;

	const int TOTAL = N(state, action);
	int n, min = __INT_MAX__, max = 0;

	for(int s = 0; s < n_states; ++s)
	{
		n = trials[s];
		if(n < min)
		{
			min = n;
		}
		if(n > max)
		{
			max = n;
		}
	}
	LAMBDA(state, action) = (max + lambdaParameter) / (TOTAL * (1 + lambdaParameter))
		- (min + lambdaParameter) / (TOTAL * (1 + lambdaParameter));


	// o Update

	for(int s = 0; s < n_states; ++s)
	{
		for(int a = 0; a < n_actions; ++a)
		{
			// - Eligibility traces
			if(s == state)
			{
				if(a == action)
				{
					#ifdef REPLACING_TRACES
					E(s, a) = 1;
					#else
					E(s, a) += 1;
					#endif
				} else
				{
					#ifdef REPLACING_TRACES
					E(s, a) = 0;
					#endif
				}
			} else
			{
				#ifdef REPLACING_TRACES
				E(s, a) = gamma * LAMBDA(s, a) * E(s, a);
				#endif
			}

			// - Action-value function
			Q(s, a) = Q(s, a) + alpha * E(s, a) * delta;

			#ifndef REPLACING_TRACES
			E(s, a) = gamma * LAMBDA(s, a) * E(s, a);
			#endif
		}
	}

	return Q(next_state, next_action);
}

/// Partial TD-MB observation (can be used with eligibility traces)

real TdMb::Observe(real reward, int next_state, int next_action)
{
	if(state >= 0)
	{
		Observe(state, action, reward, next_state, next_action);
	}

	state = next_state;
	action = next_action;

	return Q(next_state, next_action);
}

/// Get an action using the current exploration policy.
/// It calls Observe as a side-effect.

int TdMb::Act(real reward, int next_state)
{
	exploration_policy->Observe(reward, next_state);
	int next_action = exploration_policy->SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}
