/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
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

#include "FactoredPredictor.h"
#include "FactoredMarkovChain.h"
#include "BayesianPredictiveStateRepresentation.h"
#include "BayesianPredictiveStateRepresentationCTW.h"
#include <cstdlib>
#include <cstdio>

class Corridor
{
protected:
	int n_states;
	int observation;  ///< current observation
	int state; ///< current state
public:
	Corridor(int n_states_) : n_states(n_states_)
	{
		assert(n_states > 0);
	}
	void Reset()
	{
		state = 0;
		observation = 0;
	}
	int getObservation()
	{
		return observation;
	}
	bool Act(int action)
	{
		switch(action) {
		case 0:
			state--; 
			break;
		case 1:
			state++;
			break;
		default:
			break;
		}
		observation = 0;
		if (state < 0) {
			state = 0;
			observation = 1;
		} else if (state >= n_states) {
			state = n_states - 1;
			observation = 1;
		}
		return 0;
	}
};

int main(int argc, char** argv)
{

	int n_actions = 2;
	int n_obs = 2;
	
	if (argc != 4) {
		fprintf (stderr, "Usage: factored_model T states depth\n");
		return -argc;
	}
	
	int T = atoi(argv[1]);
	int n_states = atoi(argv[2]);
	int max_depth = atoi(argv[3]);

	FactoredPredictor* factored_predictor; 
	//factored_predictor = new FactoredMarkovChain(n_actions, n_obs, max_depth);

	//factored_predictor = new BayesianPredictiveStateRepresentation(n_obs, n_actions,  max_depth, 0.5);
	factored_predictor = new BayesianPredictiveStateRepresentationCTW(n_obs, n_actions,  max_depth, 0.5);
	
	Corridor environment(n_states);

	//printf ("short: %d int: %d long: %d long int: %d long long: %d\n", sizeof(short), sizeof(int), sizeof(long), sizeof (long int), sizeof (long long));
	environment.Reset();
	factored_predictor->Observe(environment.getObservation());
	for (int t=0; t<T; ++t) {
		int action = rand()%2;
		environment.Act(action);
		int observation = environment.getObservation();
		real p = factored_predictor->Observe(action, observation);
		printf ("%f\n", p);
	}


	return 0;
}
#endif
