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

#define REPLACING_TRACES

#include "TdTest.h"
#include "SmartAssert.h"

TdTest::TdTest(const int n_states_,
	const int n_actions_,
	const real gamma_,
	const real lambda_,
	const real alpha_,
	VFExplorationPolicy * const exploration_policy_,
	const real initial_value_,
	const real baseline_)
	: N_STATES(n_states_),
	N_ACTIONS(n_actions_),
	ALPHA(alpha_),
	GAMMA(gamma_),
	DEFAULT_LAMBDA(lambda_),
	exploration_policy(exploration_policy_),
	INITIAL_VALUE(initial_value_),
	BASELINE(baseline_),
	Q(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	E(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	LAMBDA(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	N(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	U(0.5),
	MEANS(n_states_, n_actions_, Matrix::CHECK_BOUNDS)
{
	make_assert(ALPHA >= 0 && ALPHA <= 1);
	make_assert(GAMMA >= 0 && GAMMA <= 1);
	make_assert(DEFAULT_LAMBDA >= 0 && DEFAULT_LAMBDA <= 1);

	Reset();

	exploration_policy->setValueMatrix(&Q);
}

TdTest::~TdTest() { }

void TdTest::Reset()
{
	state = -1;
	action = -1;

	//episode();

	for(int s = 0; s < N_STATES; ++s)
	{
		for(int a = 0; a < N_ACTIONS; ++a)
		{
			Q(s, a) = INITIAL_VALUE;

			E(s, a) = 0.0;
			LAMBDA(s, a) = DEFAULT_LAMBDA;

			vector<int> sa_trials(N_STATES, 0);
			pair<int, int> sa(s, a);
			TRIALS[sa] = sa_trials;

			N(s, a) = 0;
			MEANS(s, a) = 0;
		}
	}

}

/// Full TD-MB observation

real TdTest::Observe(const int state, const int action, const real reward, const int next_state, const int next_action)
{
	// o Update the statistics

	const pair<int, int> SA(state, action);
	vector<int> trials = TRIALS[SA];

	trials[next_state] += 1;
	N(state, action) += 1;


	// o Compututation

	// - Probabilities to go to the next state s' knowing state s and action a
	const int TOTAL = N(state, action);
	int nTrials, min = __INT_MAX__, max = 0;

	// - Means and values
	vector<double> p(N_STATES, 0), V(N_STATES, 0);
	double mean = 0;

	for(int s = 0; s < N_STATES; ++s)
	{
		// Find the min and max numbers of trials
		nTrials = trials[s];
		if(nTrials < min)
		{
			min = nTrials;
		}
		if(nTrials > max)
		{
			max = nTrials;
		}

		// Get the probability that the next state is s (from pair SA)
		p[s] = GetProbability(nTrials, TOTAL);

		// Get the value of the state s
		V[s] = QMax(s);

		// Compute the mean value of the pair SA
		mean += p[s] * V[s];
	}

	MEANS(state, action) = mean;

	// Get the variance and the standard deviation of value of the pair SA
	double variance = 0;
	for(int s = 0; s < N_STATES; ++s)
	{
		variance += p[s] * pow(V[s] - mean, 2);
	}
	double deviation = sqrt(variance);


	// o Updates

	// - Lambda (assuming positive rewards)
	double k = ((double) (max - min)) / ((double) (TOTAL * (1 + U)));
	LAMBDA(state, action) = (mean == 0) ? DEFAULT_LAMBDA : ((deviation > mean) ? k : 1 - deviation / mean);

	// - Delta
	//const real DELTA = (reward - BASELINE) + GAMMA * Q(next_state, next_action) - Q(state, action);
	const real DELTA = (reward - BASELINE) + GAMMA * MEANS(next_state, next_action) - mean;

	// - Action-value function Q for every state s and action a
	for(int s = 0; s < N_STATES; ++s)
	{
		for(int a = 0; a < N_ACTIONS; ++a)
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
				E(s, a) *= GAMMA * LAMBDA(s, a);
				#endif
			}

			// - Action-value function
			Q(s, a) += ALPHA * E(s, a) * DELTA;

			#ifndef REPLACING_TRACES
			E(s, a) *= GAMMA * LAMBDA(s, a);
			#endif
		}
	}

	return Q(next_state, next_action);
}

/// Partial TD-MB observation (can be used with eligibility traces)

real TdTest::Observe(const real reward, const int next_state, const int next_action)
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

int TdTest::Act(const real reward, const int next_state)
{
	exploration_policy->Observe(reward, next_state);
	const int next_action = exploration_policy->SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}
