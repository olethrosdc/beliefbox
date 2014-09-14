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

#ifndef UCT_MONTE_CARLO_TREE_SEARCH_LEVEL_H
#define UCT_MONTE_CARLO_TREE_SEARCH_LEVEL_H

#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
//#include "RandomNumberFile.h"

#include "real.h"
#include "Matrix.h"
#include "Grid.h"
#include "Environment.h"
#include "Random.h"
#include <limits>

template <class S, class A>
    class UCTMCL
{
 private:
  real gamma;                    // Discount factor.
  real c_uct;
  ContinuousStateEnvironment* environment; // The environment model.
  RandomNumberGenerator* rng;    // Random Number generator.
  EvenGrid discretize;
  real learning_rate; // Learning rate.
  int nBins;  // Number of grids per state space dimension.
  int MaxDepth;                  // Maximum tree depth.
  int NRollouts; // Number of sampled rollouts.
  int nActions;  // Number of available actions.
  std::vector<Matrix> Q;    // State-Actions values.
  std::vector<Matrix> C;    // Counters.
 public:
 UCTMCL(const real& gamma_, const real& c_uct_, ContinuousStateEnvironment* environment_, RandomNumberGenerator* rng_, const EvenGrid &discretize_, const real& learning_rate_, const int& MaxDepth_, const int& NRollouts_)
   : gamma(gamma_),
    c_uct(c_uct_),
    environment(environment_),
    rng(rng_),
    discretize(discretize_),
    learning_rate(learning_rate_),
    MaxDepth(MaxDepth_),
    NRollouts(NRollouts_)
    {
      nActions = environment->getNActions();    
      printf("Constructor\n");
      Q.resize(MaxDepth);
      C.resize(MaxDepth);
      for(int i = 0; i<MaxDepth; ++i) {
	Q[i] = Matrix(discretize.getNIntervals(), nActions);
	C[i] = Matrix(discretize.getNIntervals(), nActions);
      }
    };
  ~UCTMCL() {}
  // Reset statistics
  void UCT_Reset() {
    for(int i = 0; i<MaxDepth; ++i) {
      Q[i] = Matrix(discretize.getNIntervals(), nActions);
      C[i] = Matrix(discretize.getNIntervals(), nActions);
    }
  };

  double rollOut(S state_) {
    int t = 0;
    int horizon = 1000; // tree.MaxDepth - depth;
    environment->Reset();
    environment->setState(state_);

    bool running = true;
    real discount = 1.0;
    real total_reward = 0.0;
    real discounted_reward = 0.0;
    
    real reward = 0.0;
    do {
      //get current state
      S state = environment->getState();
              
      //choose an action using Random policy
      A action = rng->discrete_uniform(nActions);

      //execute the selected action
      running = environment->Act(action);
              
      // get reward
      reward = environment->getReward();
              
      total_reward += reward;
      discounted_reward += discount * reward;
              
      discount *= gamma;

      ++t;              
      if (t >= horizon) {
	running = false;
      }
    }while(running);
    //printf("Gamma = %f, discounted_reward = %f\n",tree.gamma,discounted_reward);
    return discounted_reward;
  }

  real UCT_Search(S state, int depth, bool terminal){

    
    if( terminal) {
      return 0.0;
    } else if(depth >= MaxDepth)  {
      return rollOut(state);
    }
    int index = discretize.getInterval(state);
    real SampleReturn = 0;

    int bestAction = 0;
    // printf("bestAction index = %d, value = %f \n", index,Q[depth](index, bestAction));
    real bestValue = Q[depth](index, bestAction) + 2*c_uct*sqrt( log(C[depth].RowSum(index)) / (C[depth](index,bestAction)+1) );
    
    for(int a = 1; a < nActions; ++a) {
      real value = Q[depth](index, a) + 2*c_uct*sqrt( log(C[depth].RowSum(index)) / (C[depth](index,a)+1) );
      if(value > bestValue) {
    	bestValue = value;
    	bestAction = a;
      }
    }
    environment->Reset();

    bool running = environment->Act(bestAction);
    S nextState  = environment->getState();

    real reward  = environment->getReward();
   
    SampleReturn = reward + gamma*UCT_Search(nextState, depth + 1, !running);
    C[depth](index, bestAction) += 1;
    Q[depth](index, bestAction) += (SampleReturn - Q[depth](index, bestAction))/ C[depth](index, bestAction);
    // Q(index, bestAction) = (1-learning_rate)*Q(index, bestAction) + learning_rate*SampleReturn;

    return SampleReturn;
  };

  int PlanPolicy(const S& state) {
    int rollouts = 0;
    //   printf("PlanPolicy\n");
    do {
      UCT_Search(state, 0, false);
      environment->Reset();
      environment->setState(state);
      rollouts++;
    } while(rollouts < NRollouts);

    int bestAction = 0;
    int index = discretize.getInterval(state);
    real bestValue = Q[0](index,0);
    printf("Q value = %f\n",Q[0](index,0));
    for(int a = 1; a < nActions; ++a) {
      real value = Q[0](index,a);
      printf("Q value = %f\n",Q[0](index,a));
      if(value > bestValue) {
	bestValue = value;
	bestAction = a;
      }
    }
    return bestAction;
  };
};

#endif
