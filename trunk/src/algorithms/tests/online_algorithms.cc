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
#include "PolicyEvaluation.h"
#include "ValueIteration.h"
#include "RandomMDP.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"
#include "Environment.h"
#include "ExplorationPolicy.h"
#include "Sarsa.h"

real EvaluateAlgorithm(int n_iterations,
                       OnlineAlgorithm<int,int>* algorithm,
                       DiscreteEnvironment* environment);

int main (void)
{
    int n_actions = 4;
    int n_states = 16;
    real gamma = 0.99;
    real lambda = 0.9;
    real alpha = 0.1;
    real randomness = 0.0;
    real pit_value = -100.0;
    real goal_value = 0.0;
    real step_value = -0.1;
    real epsilon = 0.1;
    int n_iterations;

    ExplorationPolicy* exploration_policy = NULL;



    

    
    exploration_policy = new EpsilonGreedy(n_actions, epsilon);

    
    OnlineAlgorithm<int, int>* algorithm;
    algorithm = new Sarsa(n_states,
                          n_actions,
                          gamma,
                          lambda,
                          alpha,
                          exploration_policy);

    DiscreteEnvironment* environment;
    environment = new RandomMDP (n_actions,
                                 n_states,
                                 randomness,
                                 step_value,
                                 pit_value,
                                 goal_value);

    //const DiscreteMDP* mdp = environment->getMDP();
    //assert(n_states == mdp->GetNStates());
    //assert(n_actions == mdp->GetNActions());

    EvaluateAlgorithm(n_iterations, algorithm, environment);
    return 0;
}

real EvaluateAlgorithm(int n_iterations,
                       OnlineAlgorithm<int, int>* algorithm,
                       DiscreteEnvironment* environment)
{
 
   environment->Reset();
    for (int iter=0; iter < n_iterations; ++iter) {
        int state = environment->getState();
        real reward = environment->getReward();
        int action = algorithm->Act(reward, state);
        bool action_ok = environment->Act(action);
        if (!action_ok) {
            environment->Reset();
        }
        std::cout << "# iter: " << iter << std::endl;
    }
    return 0.0;
}

#endif
