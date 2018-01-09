/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "MountainCar.h"
#include "CartPole.h"
#include "Pendulum.h"
#include "MersenneTwister.h"
#include "DiscreteChain.h"
#include "EasyClock.h"
#include "BasisSet.h"
#include "RepresentativeStateValueIteration.h"

int main (int argc, char** argv)
{
    MersenneTwisterRNG rng;
    //CartPole environment;
    //Pendulum environment;
    MountainCar environment;
	real gamma = 0.95;
	int grid_size = 4;
	real grid_scale = 1.0;
	int n_samples = 10;

	logmsg("Generating grid\n");
	EvenGrid grid(environment.StateLowerBound(),
				  environment.StateUpperBound(),
				  grid_size);

	logmsg("Creating kernel\n");
	// use an RBF basis for the kernel fromthe grid
	RBFBasisSet kernel(grid, grid_scale);

	logmsg("Selecting representative states\n");
	// create the set of representative states (identical to the grid)
	std::vector<Vector> states;
	std::vector<int> actions;
	for (int i=0; i<grid.getNIntervals(); ++i) {
		states.push_back(grid.getCenter(i));
	}
	// just use all the actions since it's a discrete space
	for (int i=0; i<environment.getNActions(); ++i) {
		actions.push_back(i);
	}

	logmsg("setting up RSVI\n");
	RepresentativeStateValueIteration<Vector, int, RBFBasisSet, Environment<Vector, int> > rsvi(gamma, kernel, states, actions, environment, n_samples);

	logmsg("Calculating approximate value function\n");
	rsvi.CalculateValues(0, 100);
	
    printf("\nDone\n");
    return 0.0;
}


#endif
