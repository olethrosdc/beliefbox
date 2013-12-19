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

#include "TransitionDistribution.h"
#include "Random.h"

void DiscreteTransitionDistribution::SetTransition(int state,
													int action,
													int next_state,
													real probability)
{	
	assert(probability >= 0 && probability <= 1);
	P[DiscreteTransition(state, action, next_state)] = probability;
}

int DiscreteTransitionDistribution::generate(int state, int action)
{
	real X = urandom();
	real sum = 0.0;
	for (int i=0; i<n_states; ++i) {
		real probability = P[DiscreteTransition(state, action, i)];
		sum += probability;
		if (X <= sum) {
			return i;
		}
	}
	Swarning("This statement should never be reached\n");
	return urandom(0, n_states);
}

real DiscreteTransitionDistribution::pdf(int state, int action, int next_state)
{
	return P[DiscreteTransition(state, action, next_state)];
}

