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

real TestRandomMDP(int n_states, int n_samples);


int main(int argc, char** argv)
{

    real discrete_mdp_error = TestRandomMDP(128, 32);
    if (discrete_mdp_error > 128) {
        
    }
    return 0;
}

real TestRandomMDP(int n_states, int n_samples)
{
	uint n_actions = 4;
	real randomness = 0.01;
	real step_value = -0.1;
	real pit_value = -10;
	real goal_value = 1.0;
	real discount_factor = 0.99;
    real accuracy = 1e-6;
	MersenneTwisterRNG rng;
    RandomMDP random_mdp(n_actions,
                         n_states,
                         randomness,
                         step_value,
                         pit_value,
                         goal_value,
                         &rng);
    
    DiscreteEnvironment& environment = random_mdp;
    DiscreteMDP* mdp = random_mdp.getMDP();
    RepresentativeStateModel<DiscreteEnvironment, int, int> representative_model(environment, n_samples, n_actions);
    
    ValueIteration VI(mdp, discount_factor, accuracy);
    representative_model.ComputeStateValues(discount_factor, accuracy);

    real total_error = 0;
    for (int i=0; i<n_states; ++i) {
        real V = VI.getValue(i);
        real V_approx = representative_model.getValue(i);
        printf ("%f %f\n", V, V_approx);
        total_error += abs(V - V_approx);
    }

    delete mdp;
    return total_error;
}


#endif
