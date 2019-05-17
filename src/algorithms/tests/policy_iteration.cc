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
#include "RiverSwim.h"

int main (void)
{
    int period = 30;
    int max_items = 10;
    real gamma = 0.99;
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
	RiverSwim river_swim;
    //RandomMDP random_mdp(1, 4, 0.01, 0.1, 0, 1, false);
    //random_mdp.AperiodicityTransform(0.99);
    //const DiscreteMDP* mdp = random_mdp.getMDP();
    //const DiscreteMDP* mdp = inventory_management.getMDP();
    //const DiscreteMDP* mdp = grid_world.getMDP();
    //const DiscreteMDP* mdp = optimistic_task.getMDP();
    //const DiscreteMDP* mdp = chain_task.getMDP();
	const DiscreteMDP* mdp = river_swim.getMDP();

    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();
    printf ("# states: %d, actions: %d\n", n_states, n_actions);
    AveragePolicyEvaluation average_policy_evaluation(NULL, mdp);
    //PolicyIteration policy_iteration(&average_policy_evaluation, mdp, gamma);
	PolicyIteration policy_iteration(mdp, gamma);
    ValueIteration value_iteration(mdp, gamma);
    //AverageValueIteration value_iteration(mdp);

    std::vector<real> Q(n_actions);

	real accuracy = 0.0;

    printf ("Policy iteration\n");
    policy_iteration.ComputeStateValues(accuracy, 1000);
    for (int s=0; s<mdp->getNStates(); s++) {
        printf ("%.2f ", policy_iteration.getValue(s));
    }
    printf(", D = %.1f, b = %.1f\n",
		   policy_iteration.Delta,
		   policy_iteration.baseline);
	policy_iteration.policy->Show();
	for (int s=0; s<mdp->getNStates(); s++) {
		printf("%f %f %f\n",
			   policy_iteration.getValue(s),
			   policy_iteration.getValue(s, 0),
			   policy_iteration.getValue(s, 1));
	}
	
    printf ("Value iteration\n");
    //value_iteration.ComputeStateActionValues(accuracy, 1000);
    value_iteration.ComputeStateValues(accuracy, 100000);
    for (int s=0; s<mdp->getNStates(); s++) {
        printf ("%.2f ", value_iteration.getValue(s));
    }
    printf(", D = %.1f, b = %.1f\n",
		   value_iteration.Delta,
		   value_iteration.baseline);
	FixedDiscretePolicy* value_iteration_policy = value_iteration.getPolicy();
	value_iteration_policy->Show();
	value_iteration.getStateValues().print(stdout);
	value_iteration.getValues().print(stdout);
	delete value_iteration_policy;
	{
		printf("Policy Evaluation of VI\n");
		FixedDiscretePolicy* vi_policy = value_iteration.getPolicy();
		PolicyEvaluation policy_evaluation (vi_policy, mdp, gamma);
		policy_evaluation.ComputeStateValues(accuracy, 10000);
		for (int s=0; s<mdp->getNStates(); s++) {
			printf ("%f ", policy_evaluation.getValue(s));
		}
		printf("\n");
		policy_evaluation.ComputeStateValuesFeatureExpectation(accuracy, 10000);
		for (int s=0; s<mdp->getNStates(); s++) {
			printf ("%f ", policy_evaluation.getValue(s));
		}
		printf("\n");
	}
	{
		printf("Policy Evaluation of PI\n");
		PolicyEvaluation policy_evaluation (policy_iteration.policy, mdp, gamma);
		policy_evaluation.ComputeStateValues(accuracy, 10000);
		for (int s=0; s<mdp->getNStates(); s++) {
			printf ("%f ", policy_evaluation.getValue(s));
		}
		printf("\n");
		policy_evaluation.ComputeStateValuesFeatureExpectation(accuracy, 10000);
		for (int s=0; s<mdp->getNStates(); s++) {
			printf ("%f ", policy_evaluation.getValue(s));
		}
		printf("\n");
	}
	
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
