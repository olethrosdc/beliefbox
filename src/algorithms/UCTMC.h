/* -*- Mode: C++; -*- */
// copyright (c) 2014 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef UCT_MONTE_CARLO_TREE_SEARCH_H
#define UCT_MONTE_CARLO_TREE_SEARCH_H

#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFile.h"

#include "real.h"
#include "Matrix.h"
#include "Grid.h"
#include "Environment.h"
#include "Random.h"
#include <limits>

template <class S, class A>
class UCTMC
{
private:
	real gamma;                    // Discount factor.
	real c_uct;
	ContinuousStateEnvironment* environment; // The environment model.
	RandomNumberGenerator* rng;    // Random Number generator.
	EvenGrid discretize;
	real learning_rate; // Learning rate.
	real lambda; // mixing
	int nBins;  // Number of grids per state space dimension.
	int MaxDepth;                  // Maximum tree depth.
	int NRollouts; // Number of sampled rollouts.
	int nActions;  // Number of available actions.
	Matrix Q;    // State-Actions values.
	Matrix C;    // Counters.
public:
	UCTMC(const real& gamma_, const real& c_uct_, ContinuousStateEnvironment* environment_, RandomNumberGenerator* rng_, const EvenGrid &discretize_, const real& learning_rate_, const real& lambda_, const int& MaxDepth_, const int& NRollouts_)
		: gamma(gamma_),
		  c_uct(c_uct_),
		  environment(environment_),
		  rng(rng_),
		  discretize(discretize_),
		  learning_rate(learning_rate_),
		  lambda(lambda_),
		  MaxDepth(MaxDepth_),
		  NRollouts(NRollouts_)
	{
		nActions = environment->getNActions();    

		Q = Matrix(discretize.getNIntervals(), nActions);
		C = Matrix(discretize.getNIntervals(), nActions);
	};
	~UCTMC() {}
	// Reset statistics
	void UCT_Reset() {
		Q = Matrix(discretize.getNIntervals(), nActions);
		C = Matrix(discretize.getNIntervals(), nActions);
	};

	real UCT_Search(S state, int depth, bool terminal){

		if(depth >= MaxDepth ||  terminal) {
			return 0.0;
		}
		int index = discretize.getInterval(state);
		real SampleReturn = 0;

		int bestAction = 0;
		real bestValue = Q(index, bestAction) + 2*c_uct*sqrt( log(C.RowSum(index)) / (C(index,bestAction)+1) );
    
		for(int a = 1; a < nActions; ++a) {
			real value = Q(index, a) + 2*c_uct*sqrt( log(C.RowSum(index)) / (C(index,a)+1) );
			if(value > bestValue) {
				bestValue = value;
				bestAction = a;
			}
		}
		//environment->Reset();
    
		bool running = environment->Act(bestAction);
		S nextState  = environment->getState();

		real reward  = environment->getReward();
   
		SampleReturn = reward + gamma*UCT_Search(nextState, depth + 1, !running);
		C(index, bestAction) += 1;
		Q(index, bestAction) += (SampleReturn - Q(index, bestAction))/ C(index, bestAction);
		// Note: it's best to use learning-rate when lambda < 1.
	
		//Q(index, bestAction) = (1-learning_rate)*Q(index, bestAction) + learning_rate*SampleReturn;
		real LambdaReturn = lambda * SampleReturn + (1 - lambda) * Q(index, bestAction);
		return LambdaReturn;
	};
	int PlanPolicy(const S& state) {
		int rollouts = 0;

		do {
  
			UCT_Search(state, 0, false);
			environment->Reset();
			environment->setState(state);
			rollouts++;
		} while(rollouts < NRollouts);

		int bestAction = 0;
		int index = discretize.getInterval(state);
		real bestValue = Q(index,0);

		for(int a = 1; a < nActions; ++a) {
			real value = Q(index,a);
			if(value > bestValue) {
				bestValue = value;
				bestAction = a;
			}
		}
		return bestAction;
	};
};

#endif
