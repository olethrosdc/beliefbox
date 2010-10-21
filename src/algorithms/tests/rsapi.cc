/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "RSAPI.h"
#include "Rollout.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "RandomPolicy.h"
#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"

int main(void)
{
	MersenneTwisterRNG rng;

	// Create a new environment
	Environment<Vector, int>* environment;

	//environment = new MountainCar();
	environment = new Pendulum();
	
	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;
	
	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), &rng);
	
    RSAPI rsapi(environment, &rng);

    int n_states = 100;
    
    int state_dimension = environment->getNStates();
    Vector S_L = environment->StateLowerBound();
    Vector S_U = environment->StateUpperBound();
    
    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);

    for (int k=0; k<n_states; ++k) {
        Vector state(S_L.Size());
        for (int i=0; i<S_L.Size(); ++i) {
            state(i) = rng.uniform(S_L(i), S_U(i));
        }
        rsapi.AddState(state);
        //printf("# Adding state: "); state.print(stdout);
    }
    
    rsapi.setPolicy(policy);
    rsapi.NewRandomRollouts(100, 1000);
	delete policy;
	delete environment;
}

#endif
