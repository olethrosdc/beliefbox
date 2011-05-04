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
#include "ValueIteration.h"
#include "PolicyIteration.h"
#include "AverageValueIteration.h"
#include "AveragePolicyEvaluation.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "RandomMDP.h"
#include "OptimisticTask.h"
#include "DiscreteChain.h"
#include "DiscretePolicy.h"

int main (void)
{
    int period = 30;
    int max_items = 10;
    real gamma = 0.9;
    real demand = 0.1;
    real random = 0.0;
    real pit = -100.0;
    real goal = 0.0;
    real step = -0.1;
    real margin = 1.1;
    
    //InventoryManagement inventory_management (period, max_items, demand, margin);
    
    Gridworld grid_world("/home/olethros/projects/beliefbox/dat/maze2", random, pit, goal, step);
    OptimisticTask optimistic_task(0.1, 0.01);
    DiscreteChain chain_task(5);
    //RandomMDP random_mdp(1, 4, 0.01, 0.1, 0, 1, false);
    //random_mdp.AperiodicityTransform(0.99);
    //const DiscreteMDP* mdp = random_mdp.getMDP();
    //const DiscreteMDP* mdp = inventory_management.getMDP();
    //const DiscreteMDP* mdp = grid_world.getMDP();
    //const DiscreteMDP* mdp = optimistic_task.getMDP();
    const DiscreteMDP* mdp = chain_task.getMDP();

    int n_states = mdp->GetNStates();
    int n_actions = mdp->GetNActions();
    printf ("# states: %d, actions: %d\n", n_states, n_actions);
    AveragePolicyEvaluation policy_evaluation(NULL, mdp);
    PolicyIteration policy_iteration(&policy_evaluation, mdp, gamma);
    ValueIteration value_iteration(mdp, gamma);
    //AverageValueIteration value_iteration(mdp);

    std::vector<real> Q(n_actions);

    printf ("Policy iteration\n");
    policy_iteration.ComputeStateValues(0.01, 1000);
    for (int s=0; s<mdp->GetNStates(); s++) {
        printf ("%1.f ", policy_iteration.getValue(s));
    }
    printf(", D = %.1f, b = %.1f\n",
              policy_iteration.Delta,
              policy_iteration.baseline);

    printf ("Value iteration\n");
    value_iteration.ComputeStateActionValues(0.01, 100000);
    for (int s=0; s<mdp->GetNStates(); s++) {
        printf ("%1.f ", value_iteration.getValue(s));
    }
    printf(", D = %.1f, b = %.1f\n",
              value_iteration.Delta,
              value_iteration.baseline);

    //mdp->ShowModel();
   if (1) {
       FILE* f = fopen("test.dot", "w");
       if (f) {
           fprintf (f, "digraph MDP {\n");
           fprintf (f, "ranksep=2; rankdir=LR; \n");
           mdp->dotModel(f);
           fprintf (f, "}\n");
           fclose(f);
       }
   }

   
   

}


#endif
