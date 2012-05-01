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
    int period = 30;
    int max_items = 10;
    real gamma = 0.95;
    real demand = 0.1;
    real random = 0.1;
    real pit = -100.0;
    real goal = 1.0;
    real step = -0.1;
	real margin = 1.1;
    MersenneTwisterRNG rng;
    //InventoryManagement inventory_management (period, max_items, demand, margin);
    //const DiscreteMDP* mdp = inventory_management.getMDP();

    Gridworld grid_world("/home/olethros/projects/beliefbox/dat/maze0a", random, pit, goal, step);
    const DiscreteMDP* mdp = grid_world.getMDP();
    
    //RandomMDP random_mdp(32, 32, 0.001, 0.1, 0, 1, &rng, false);
    //const DiscreteMDP* mdp = random_mdp.getMDP();
    
    //DiscreteChain chain(5);
    //const DiscreteMDP* mdp = chain.getMDP();


    
    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();

    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesStandard(0.00, 1000);
        double end_time = GetCPU();
        printf("\nStandard: %f\n", end_time - start_time);
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

    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesAsynchronous(0.00, 1000);
        double end_time = GetCPU();
        printf("\nAsynchronous: %f\n", end_time - start_time);

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

    {
        ValueIteration value_iteration(mdp, gamma);
        double start_time = GetCPU();
        value_iteration.ComputeStateValuesElimination(0.00, 1000);
        double end_time = GetCPU();
        printf("\nElimination: %f\n", end_time - start_time);

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
    printf("\nDone\n");
    return 0.0;
}


#endif
