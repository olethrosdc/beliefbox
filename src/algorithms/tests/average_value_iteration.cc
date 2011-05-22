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
//#include "GoalStateIteration.h"
#include "AverageValueIteration.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "RandomMDP.h"
#include "DiscretePolicy.h"

struct Stats
{
    real Delta;
    bool max_iter_reached;
};

Stats Evaluate();

void TestModel(const DiscreteMDP* mdp,
               const int max_iter,
               const real threshold,
               const bool relative,
               const char* fname);

DiscreteMDP TestExample8_5_1()
{
    DiscreteMDP mdp(2, 1, NULL, NULL);
    mdp.setTransitionProbability(0, 0, 0, 0.0);
    mdp.setTransitionProbability(0, 0, 1, 1.0);
    mdp.setTransitionProbability(1, 0, 0, 1.0);
    mdp.setTransitionProbability(1, 0, 1, 0.0);
    mdp.addFixedReward(0, 0, 0);
    mdp.addFixedReward(1, 0, 0);
    return mdp;
}

// for this mdp, E[r | s1, \pi] = 10, E[r | s2, \pi] = -1
// for all policies
DiscreteMDP TestExample8_5_2()
{
    // state 2 has only one action, modelling it as two identical ones
    DiscreteMDP mdp(2, 2, NULL, NULL);
    // state 1
    // action 1 r = 5
    mdp.addFixedReward(0, 0, 5);
    mdp.addFixedReward(0, 1, 10);
    mdp.addFixedReward(1, 0, -1);
    mdp.addFixedReward(1, 1, -1);

    // s1 -(a1)-> s1, p = 0.5
    mdp.setTransitionProbability(0, 0, 0, 0.5);
    mdp.setTransitionProbability(0, 0, 1, 0.5);
    // s1 -(a1)-> s2, p = 1
    mdp.setTransitionProbability(0, 1, 0, 0.0);
    mdp.setTransitionProbability(0, 1, 1, 1.0);

    // s2 -> s2, 
    mdp.setTransitionProbability(1, 0, 0, 0.0);
    mdp.setTransitionProbability(1, 0, 1, 1.0);
    mdp.setTransitionProbability(1, 1, 0, 0.0);
    mdp.setTransitionProbability(1, 1, 1, 1.0);

    return mdp;
}

// for this mdp, E[r | s1, \pi] = 10, E[r | s2, \pi] = -1
// for all policies
DiscreteMDP TestExample9_5_1()
{
    // state 2 has only one action, modelling it as two identical ones
    DiscreteMDP mdp(3, 3, NULL, NULL);
    // state 1, a_1 0,1
    mdp.addFixedReward(0, 0, 0); 
    mdp.addFixedReward(0, 1, 2); 
    mdp.addFixedReward(0, 2, 2);
    // state 2
    mdp.addFixedReward(1, 0, 1);
    mdp.addFixedReward(1, 1, 1);
    mdp.addFixedReward(1, 2, 3);
    // state 3
    mdp.addFixedReward(2, 0, 2);
    mdp.addFixedReward(2, 1, 4);
    mdp.addFixedReward(2, 2, 4);

    // a1,1 -> s2
    mdp.setTransitionProbability(0, 0, 0, 0.0);
    mdp.setTransitionProbability(0, 0, 1, 1.0);
    mdp.setTransitionProbability(0, 0, 2, 0.0);
    // a1,2 -> s1
    mdp.setTransitionProbability(0, 1, 0, 1.0);
    mdp.setTransitionProbability(0, 1, 1, 0.0);
    mdp.setTransitionProbability(0, 1, 2, 0.0);
    // a1,3 -> s1
    mdp.setTransitionProbability(0, 2, 0, 1.0);
    mdp.setTransitionProbability(0, 2, 1, 0.0);
    mdp.setTransitionProbability(0, 2, 2, 0.0);

    // a2,1 -> s3
    mdp.setTransitionProbability(1, 0, 0, 0.0);
    mdp.setTransitionProbability(1, 0, 1, 0.0);
    mdp.setTransitionProbability(1, 0, 2, 1.0);
    // a2,2 -> s1
    mdp.setTransitionProbability(1, 1, 0, 1.0);
    mdp.setTransitionProbability(1, 1, 1, 0.0);
    mdp.setTransitionProbability(1, 1, 2, 0.0);
    // a2,3 -> s2
    mdp.setTransitionProbability(1, 2, 0, 0.0);
    mdp.setTransitionProbability(1, 2, 1, 1.0);
    mdp.setTransitionProbability(1, 2, 2, 0.0);

    // a3,1 -> s2
    mdp.setTransitionProbability(2, 0, 0, 0.0);
    mdp.setTransitionProbability(2, 0, 1, 1.0);
    mdp.setTransitionProbability(2, 0, 2, 0.0);
    // a3,2 -> s3
    mdp.setTransitionProbability(2, 1, 0, 0.0);
    mdp.setTransitionProbability(2, 1, 1, 0.0);
    mdp.setTransitionProbability(2, 1, 2, 1.0);
    // a3,3 -> s3
    mdp.setTransitionProbability(2, 2, 0, 0.0);
    mdp.setTransitionProbability(2, 2, 1, 0.0);
    mdp.setTransitionProbability(2, 2, 2, 1.0);


    return mdp;
}

int main (void)
{
    int max_iter = 10;
    real threshold = 0.1;
    {
        printf("Testing example 8.5.1\n");
        DiscreteMDP mdp = TestExample8_5_1();
        TestModel(&mdp, max_iter, threshold, false, "example_8_5_1.dot");
    }
    {
        printf("\nTesting example 8.5.2\n");
        printf("---------------------\n");
        DiscreteMDP mdp = TestExample8_5_2();
        printf("\nAbsolute value iteration\n");
        TestModel(&mdp, max_iter, threshold, false, "example_8_5_2.dot");
        printf("\nRelative value iteration\n");
        TestModel(&mdp, max_iter, threshold, true, "example_8_5_2_relative.dot");
        mdp.AperiodicityTransform(0.99);
        printf("\nAbsolute value iteration\n");
        TestModel(&mdp, max_iter, threshold, false, "example_8_5_2_aperiodic_dot");
        printf("\nRelative value iteration\n");
        TestModel(&mdp, max_iter, threshold, true, "example_8_5_2_aperiodic_relative.dot");

    }
    {
        printf("\nTesting example 9.5.1\n");
        printf("---------------------\n");
        DiscreteMDP mdp = TestExample9_5_1();
        printf("\nAbsolute value iteration\n");
        TestModel(&mdp, max_iter, threshold, false, "example_9_5_1.dot");
        printf("\nRelative value iteration\n");
        TestModel(&mdp, max_iter, threshold, true, "example_9_5_1_relative.dot");
        mdp.AperiodicityTransform(0.99);
        printf("\nAbsolute value iteration\n");
        TestModel(&mdp, max_iter, threshold, false, "example_9_5_1_aperiodic_dot");
        printf("\nRelative value iteration\n");
        TestModel(&mdp, max_iter, threshold, true, "example_9_5_1_aperiodic_relative.dot");

    }

#if 0
    real mean_Delta = 0.0;
    int failures = 0;
    int iter_failures = 0;
    int n_tests = 10;
    for (int test=0; test<n_tests; test++) {
        Stats stats = Evaluate();
        mean_Delta += stats.Delta;
        if (stats.Delta > 1.0) {
            failures++;
        }
        if (stats.max_iter_reached) {
            iter_failures++;
        }
    }
    mean_Delta /= (real) n_tests;
    printf ("Average Delta: %f, Delta Failures: %.2f%%, Iter Failures: %.2f%%\n",
            mean_Delta,
            ((real) (100*failures))/((real) n_tests),
            ((real) (100*iter_failures))/((real) n_tests));
    return failures;
#else
    return 0;
#endif
}





Stats Evaluate() {
    real random = 0.0;
    real pit = -1.0;
    real goal = 1.0;
    real step = -0.1;
    //InventoryManagement inventory_management (period, max_items, demand);

    //Gridworld grid_world("maze2", width, height, 4, random, pit, goal, step);
    RandomMDP random_mdp(2, 4, random, step, pit, goal, false);
    const DiscreteMDP* mdp = random_mdp.getMDP();
    //GoalStateIteration goal_state_iteration(mdp);
    real diameter = 5.0;//goal_state_iteration.GetMDPDiameter(0.0001, 1000);
    real tau = 0.9;

    random_mdp.AperiodicityTransform(tau);
    mdp->Check();
    
    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();
    
    std::vector<Vector> p(n_states);
    for (int s=0; s<n_states; s++) {
        p[s].Resize(n_actions);
        for (int a=0; a<n_actions; a++) {
            p[s][a] = 1.0 / (real) n_actions;
            //printf ("R(%d,%d) = %f\n", s, a, mdp->getExpectedReward(s, a));
        }
    }

    AverageValueIteration value_iteration(mdp, true, true);

    std::vector<real> Q(n_actions);
    
    value_iteration.ComputeStateValues(0.1, 1000);
#if 0
    for (int s=0; s<mdp->getNStates(); s++) {
        printf ("%1.f ", value_iteration.getValue(s));
    }
    printf("\n");
#endif
    
    if (value_iteration.Delta > 100) {// || value_iteration.max_iter_reached) {
        printf ("D=%f, diameter = %f", value_iteration.Delta, diameter);
        if (value_iteration.max_iter_reached == true) {
            printf(" TERM\n");
        } else {
            printf("\n");
        }
        FILE* f = fopen("test.dot", "w");
        if (f) {
            fprintf (f, "digraph MDP {\n");
            fprintf (f, "ranksep=2; rankdir=LR; \n");
            mdp->dotModel(f);
            for (int s=0; s<mdp->getNStates(); s++) {
                fprintf (f,
                         "s%d [label = \"%.2f\"];\n",
                         s,
                         value_iteration.getValue(s));
            }
            fprintf (f, "}\n");
            fclose(f);
        }
    }

    Stats stats = {value_iteration.Delta, value_iteration.max_iter_reached};
    return stats;
}

void TestModel(const DiscreteMDP* mdp,
               const int max_iter,
               const real threshold,
               const bool relative,
               const char* fname)
{
    const int n_states = mdp->getNStates();
    const int n_actions = mdp->getNActions();
    mdp->Check();

    // Get MDP diameter
    //GoalStateIteration goal_state_iteration(mdp);
    real diameter = 5.0;//goal_state_iteration.GetMDPDiameter(threshold, max_iter);
    
    
    // Do value iteration
    AverageValueIteration value_iteration(mdp, relative, true);

    printf ("State values\n");
    for (int i=0; i<max_iter; i++) {
        value_iteration.ComputeStateValues(threshold, 1);
        printf ("%d : ", i);
        for (int s=0; s<mdp->getNStates(); s++) {
            printf ("%+4.4f ", value_iteration.getValue(s));
        }
        printf("%f\n", value_iteration.Delta);
        //if (i > 0 && value_iteration.Delta < threshold) {
        //break;
        //}
    }

    
    printf ("D=%f, diameter = %f\n", value_iteration.Delta, diameter);
    
    if (fname) {
        FILE* f = fopen(fname, "w");
        if (f) {
            fprintf (f, "digraph MDP {\n");
            fprintf (f, "ranksep=2; rankdir=LR; \n");
            mdp->dotModel(f);
            for (int s=0; s<mdp->getNStates(); s++) {
                fprintf (f,
                         "s%d [label = \"%.2f\"];\n",
                         s,
                         value_iteration.getValue(s));
            }
            fprintf (f, "}\n");
            fclose(f);
        }
    }
}



#endif
