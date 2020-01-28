// -*- Mode: c++ -*-
// copyright (c) 2007-2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SimpleTransitionDistribution.h"
#include "Random.h"

SimpleTransitionDistribution::~SimpleTransitionDistribution()
{
if (arr)
	free(arr);		// << Somehow this is freeing up automatically
/*
	int n_pairs = 0;
	for (int s=0; s<n_states; s++) {	
		for (int a=0; a<n_actions; a++) {
			DiscreteStateAction SA(s, a);
			auto got = next_states.find(SA);
			if (got != next_states.end()) {
				n_pairs++;
			}
		}
	}
*/
	//printf("%d transitions, %d state-actions, %d state-action pairs saved\n",
	//	   P.size(),
	//next_states.size(),
	//n_pairs);
}

void SimpleTransitionDistribution::SetTransition(int state,
													int action,
													int next_state,
													real probability)
{	
	assert(probability >= 0 && probability <= 1);
//	DiscreteTransition transition = DiscreteTransition(state, action, next_state);
	if (probability > 0) {
//		P[transition] = probability;						<<< As I guessed, just commenting this out works superfast
		arr[(state+action*n_states)*n_states + next_state] = probability;
		DiscreteStateAction SA(state, action);
		next_states[SA].insert(next_state);
	} else {
		// erase transition
//		P.erase(transition);
		arr[(state+action*n_states)*n_states + next_state] = 0.0f;
		// erase state from the set of next states
		DiscreteStateAction SA(state, action);
		auto got = next_states.find(SA);
		if (got != next_states.end()) {
			DiscreteStateSet states = got->second;
			states.erase(next_state);
		}
	}
}

real SimpleTransitionDistribution::GetTransition(int state,
												   int action,
												   int next_state) const
{	
	return arr[(state+action*n_states)*n_states + next_state];
}

int SimpleTransitionDistribution::generate(int state, int action) const
{
	real X = urandom();
	real sum = 0.0;
	for (int i=0; i<n_states; ++i) {
		real probability = GetTransition(state, action, i);
		sum += probability;
		if (X <= sum) {
			return i;
		}
	}
	Swarning("This statement should never be reached\n");
	return urandom(0, n_states);
}

real SimpleTransitionDistribution::pdf(int state, int action, int next_state) const
{
	return GetTransition(state, action, next_state);
}

