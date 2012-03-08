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
#include "NormalDistribution.h"

TdBma::TdBma(int n_states_,
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
	N(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Means(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Variances(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
	Deviations(n_states_, n_actions_, Matrix::CHECK_BOUNDS)
{
	assert(lambda >= 0 && lambda <= 1);
	assert(alpha >= 0 && alpha <= 1);
	assert(gamma >= 0 && gamma <= 1);

	Reset();

	exploration_policy->setValueMatrix(&Q);
}

TdBma::~TdBma() { }

void TdBma::Reset()
{
	state = -1;
	action = -1;

	T = 20;
	t = 0;

	//episode();

	for(int s = 0; s < n_states; s++)
	{
		for(int a = 0; a < n_actions; a++)
		{
			Q(s, a) = initial_value;
			N(s, a) = 0.0;
			Means(s, a) = 0.0;
			Variances(s, a) = 0.0;
		}
	}

}

real TdBma::Observe(int state, int action, real reward, int next_state, int next_action)
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
		real rBma, nrBma = nr;

		// - Q-values for each method (TD and MC)
		real qTd, qMc;

		// - The normal distribution N(μ, σ^2)
		//   and intermediary variables to compute lambda
		real nTd, nMc;

		// - Intermediary variable to update μ and σ
		real delta, meanSquared;


		// o Browse the episode
		while(!episode.empty())
		{
			// o Get the state, action and reward of this step
			pair< pair< int, int >, real > sar = episode.back();
			s = sar.first.first;
			a = sar.first.second;
			r = sar.second;
			//printf("State = %d, Action = %d => Reward = %lf\n", s, a, r);


			// o Compute the model expectations
			if(t < T)
			{
				// - Compute the Q-values for TD and MC methods
				qTd = nr + gamma * Q(ns, na);
				qMc = nr + gamma * nrBma;

				// - Compute lambda
				// The normal distribution N(μ, σ^2)
				NormalDistribution nd(Means(s, a), Deviations(s, a));
				nMc = pow(nd.pdf(qMc), N(s, a));
				nTd = pow(nd.pdf(qTd), N(s, a));
				lambda = nMc / (nMc + nTd);

				// - Compute the BMA return
				rBma = (1 - lambda) * qTd + lambda * qMc;
			}


			// o Update the action-value function Q
			Q(s, a) = Q(s, a) + alpha * (rBma - Q(s, a));


			// o Update μ, σ sufficient statistics (Knuth, 1998)
			N(s, a) += 1;
			delta = r - Means(s, a);
			Means(s, a) += delta / N(s, a);
			meanSquared = delta * (r - Means(s, a));
			if(N(s, a) > 1)
			{
				Variances(s, a) += meanSquared / N(s, a);
				Deviations(s, a) = sqrt(Variances(s, a));
			} else
			{
				Variances(s, a) += meanSquared;
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

real TdBma::Observe(real reward, int next_state, int next_action)
{
	if(state >= 0)
	{
		Observe(state, action, reward, next_state, next_action);
	}

	state = next_state;
	action = next_action;

	return Q(next_state, next_action);
}

int TdBma::Act(real reward, int next_state)
{
	exploration_policy->Observe(reward, next_state);
	int next_action = exploration_policy->SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}
