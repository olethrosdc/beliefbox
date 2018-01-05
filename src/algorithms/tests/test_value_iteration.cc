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
#include "PolicyGradient.h"
#include "AverageValueIteration.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "RandomMDP.h"
#include "DiscretePolicy.h"
#include "MersenneTwister.h"
#include "DiscreteChain.h"
#include "EasyClock.h"

int main (void)
{
    MersenneTwisterRNG rng;
#if 0
    int period = 30;
    int max_items = 10;
    real demand = 0.1;
    real random = 0.1;
    real pit = -100.0;
    real goal = 1.0;
    real step = -0.1;
    real margin = 1.1;

    InventoryManagement inventory_management (period, max_items, demand, margin);
    Gridworld grid_world("/home/olethros/projects/beliefbox/dat/maze9", random, pit, goal, step);
    RandomMDP random_mdp(32, 32, 0.001, 0.1, 0, 1, &rng, false);


#endif
    
    DiscreteChain chain(5);
    const DiscreteMDP* mdp = chain.getMDP();

    //const DiscreteMDP* mdp = inventory_management.getMDP();
    //const DiscreteMDP* mdp = grid_world.getMDP();
    //const DiscreteMDP* mdp = random_mdp.getMDP();
    

    real gamma = 0.95;    
    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();
    int n_iterations = 1000;
    real accuracy = 1e-9;

    bool test_synchronous = true;
    bool test_asynchronous = true;
    bool test_elimination = true;
    bool test_gradient = true;
    
    if (test_synchronous)
    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesStandard(accuracy, n_iterations);
        double end_time = GetCPU();
        printf("\nStandard time: %f\n", end_time - start_time);
        FixedDiscretePolicy* policy = value_iteration.getPolicy();
        for (int s=0; s<n_states; ++s) {
            printf (" %d ", ArgMax(policy->getActionProbabilitiesPtr(s)));
        }
        printf("\n");
        for (int s=0; s<n_states; ++s) {
            printf (" %.1f ", value_iteration.getValue(s));
        }
        delete policy;
    }

    if (test_asynchronous)
    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesAsynchronous(accuracy, n_iterations);
        double end_time = GetCPU();
        printf("\nAsynchronous time: %f\n", end_time - start_time);

        FixedDiscretePolicy* policy = value_iteration.getPolicy();
        for (int s=0; s<n_states; ++s) {
            printf (" %d ", ArgMax(policy->getActionProbabilitiesPtr(s)));
        }
        printf("\n");
        for (int s=0; s<n_states; ++s) {
            printf (" %.1f ", value_iteration.getValue(s));
        }
        delete policy;
    }

    if (test_elimination)
    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesElimination(accuracy, n_iterations);
        double end_time = GetCPU();
        printf("\nElimination time: %f\n", end_time - start_time);

        FixedDiscretePolicy* policy = value_iteration.getPolicy();
        for (int s=0; s<n_states; ++s) {
            printf (" %d ", ArgMax(policy->getActionProbabilitiesPtr(s)));
        }
        printf("\n");
        for (int s=0; s<n_states; ++s) {
            printf (" %.1f ", value_iteration.getValue(s));
        }
        delete policy;
    }

    
    if (test_gradient)
    {
        real step_size = 0.1;
        PolicyGradient policy_gradient(mdp, gamma, step_size);
        double start_time = GetCPU();
        policy_gradient.ModelBasedGradient(accuracy, n_iterations);
        double end_time = GetCPU();
        printf("\nGradient time time: %f\n", end_time - start_time);

        FixedDiscretePolicy* policy = policy_gradient.getPolicy();
        for (int s=0; s<n_states; ++s) {
            Vector* pi = policy->getActionProbabilitiesPtr(s);
            int a = ArgMax(pi);
            printf (" %d@%.1f ", a, (*pi)(a));
        }
        printf("\n");
        for (int s=0; s<n_states; ++s) {
            printf (" %.1f ", policy_gradient.getValue(s));
        }
    }
    

    
    printf("\nDone\n");
    return 0.0;
}


#endif
