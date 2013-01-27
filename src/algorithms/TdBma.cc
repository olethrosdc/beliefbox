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

#include "TdBma.h"
#include "SmartAssert.h"

TdBma::TdBma(const int n_states_,
	const int n_actions_,
	const real gamma_,
	const real lambda_,
	const real alpha_,
	VFExplorationPolicy* exploration_policy_,
	const real initial_value_,
	const real baseline_)
	: N_STATES(n_states_),
	N_ACTIONS(n_actions_),
	ALPHA(alpha_),
	GAMMA(gamma_),
	lambda(lambda_),
	exploration_policy(exploration_policy_),
	INITIAL_VALUE(initial_value_),
	BASELINE(baseline_),
	Q(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	N(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Means(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Variances(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Deviations(n_states_, n_actions_, Matrix::CHECK_BOUNDS)
{
	assert(ALPHA >= 0 && ALPHA <= 1);
	assert(GAMMA >= 0 && GAMMA <= 1);
	assert(lambda >= 0 && lambda <= 1);

	Reset();

	exploration_policy->setValueMatrix(&Q);
}

TdBma::~TdBma() { }

void TdBma::Reset()
{
	state = -1;
	action = -1;

	T = (int) (log(N_STATES) / log(N_ACTIONS));
	t = 0;

	//episode();

	for(int s = 0; s < N_STATES; s++)
	{
		for(int a = 0; a < N_ACTIONS; a++)
		{
			Q(s, a) = INITIAL_VALUE;
			N(s, a) = 0.0;
			Means(s, a) = 0.0;
			Variances(s, a) = 0.0;
		}
	}

}

/// Full TD-MB observation (no eligibility traces)

real TdBma::Observe(const int state, const int action, const real reward, const int next_state, const int next_action)
{
	// o If the episode is finished (episode = T steps),
	//   use it to refine the estimate Q
	if(t == T)
	{
		// o Declare and init variables
		// - State, action, next state and next action
		int s, a, ns = state, na = action;

		// - Reward and next reward
		real r, nr = reward;

		// - BMA return and next BMA return
		real rBma = nr, nrBma = rBma;

		// - Q-values for each method (TD and MC)
		real qTd, qMc;

		// - The normal distribution N(μ, σ^2)
		//   and intermediary variables to compute lambda
		real nTd, nMc;

		// - Intermediary variable to update μ and σ
		real delta, meanSquared = 0;


		// o Browse the episode
		while(!episode.empty())
		{
			// o Get the state, action and reward of this step
			const pair< pair< int, int >, real > sar = episode.back();
			s = sar.first.first;
			a = sar.first.second;
			r = sar.second;
			//printf("State = %d, Action = %d => Reward = %lf\n", s, a, r);


			// o Compute the model expectations
			if(t < T)
			{
				// - Compute the Q-values for TD and MC methods
				qTd = nr + GAMMA * Q(ns, na);
				qMc = nr + GAMMA * nrBma;

				// - Compute lambda
				// The normal distribution N(μ, σ^2)
				const NormalDistribution nd(Means(s, a), Deviations(s, a));
				nMc = pow(nd.pdf(qMc), N(s, a));
				nTd = pow(nd.pdf(qTd), N(s, a));
				lambda = nMc / (nMc + nTd);

				// - Compute the BMA return
				rBma = (1 - lambda) * qTd + lambda * qMc;
			}


			// o Update the action-value function Q
			Q(s, a) += ALPHA * (rBma - Q(s, a));


			// o Update μ, σ sufficient statistics (Knuth, 1998)
			N(s, a) += 1;
			delta = rBma - Means(s, a);
			Means(s, a) += delta / N(s, a);
			meanSquared += delta * (rBma - Means(s, a));
			if(N(s, a) > 1)
			{
				Variances(s, a) = meanSquared / N(s, a);
				Deviations(s, a) = sqrt(Variances(s, a));
			}


			// o Update the next state, action, reward and BMA return
			ns = s;
			na = a;
			nr = r;
			nrBma = rBma;


			episode.pop_back();
		}

		t = 0;
	}

	// o Insert the triplet (state, action and reward) of this step
	//   into the current episode
	pair< int, int > sa(state, action);
	pair< pair< int, int >, real > sar(sa, reward);
	episode.push_back(sar);
	++t;

	return Q(next_state, next_action);
}

real TdBma::Observe(const real reward, const int next_state, const int next_action)
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

int TdBma::Act(const real reward, const int next_state)
{
	exploration_policy->Observe(reward, next_state);
	int next_action = exploration_policy->SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}
