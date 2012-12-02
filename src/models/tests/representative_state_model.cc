// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/** Test that the BPSR model predicts the next observations well.
    
    The BPSRModel gives probabilities of next observations.
    It needs a full history of observations and actions to predict
    the next observation.
    
    However, it always has observations up to x_t and actions
    up to a_{t-1}.

    The BayesianPredictiveStateRepresentation shares the same
    problem. The main difficulty is that the context in the BPSR
    is defined via Factored Markov Chains.

    The context in a FMC is the observation-action history from time
    t-D to time t. It is necessary to have the context in the FMC in
    order to find the right node in the context tree.

 */

#ifdef MAKE_MAIN

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include "Demonstrations.h"
#include "RepresentativeStateModel.h"
#include "Grid.h"
#include "RandomMDP.h"
#include "ValueIteration.h"
#include "MersenneTwister.h"
int main(int argc, char** argv)
{
	
	uint n_actions = 4;
	uint n_states = 128;
	real randomness = 0.01;
	real step_value = -0.1;
	real pit_value = -10;
	real goal_value = 1.0;
	real discount_factor = 0.99;
	MersenneTwisterRNG rng;
	RandomMDP random_mdp(n_actions,
						 n_states,
						 randomness,
						 step_value,
						 pit_value,
						 goal_value,
						 &rng);

	DiscreteEnvironment& environment = random_mdp;
	Demonstrations<int, int> demonstrations;
	FixedDiscretePolicy policy(n_states, n_actions);
	demonstrations.Simulate(environment,
							policy,
							discount_factor,
							-1);
	
}


#endif
