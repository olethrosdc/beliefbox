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

int main (void)
{
    real mean_Delta = 0.0;
    int failures = 0;
    int iter_failures = 0;
    int n_tests = 100;
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
}

Stats Evaluate() {
    int period = 30;
    int max_items = 10;
    real gamma = 1.0;
    real demand = 0.1;
    real random = 0.0;
    real pit = -100.0;
    real goal = 0.0;
    real step = -0.1;
    uint width = 8;
    uint height = 8;
    InventoryManagement inventory_management (period, max_items, demand);

    //Gridworld grid_world("maze2", width, height, 4, random, pit, goal, step);
    RandomMDP random_mdp(2, 4, 0.0, -0.1, -1, 1, false);
    const DiscreteMDP* mdp = random_mdp.getMDP();
    real tau = 0.75;
    random_mdp.AperiodicityTransform(tau);
    mdp->Check();
    //const DiscreteMDP* mdp = inventory_management.getMDP();
    //const DiscreteMDP* mdp = grid_world.getMDP();
    
    int n_states = mdp->GetNStates();
    int n_actions = mdp->GetNActions();
    
    std::vector<Vector> p(n_states);
    for (int s=0; s<n_states; s++) {
        p[s].Resize(n_actions);
        for (int a=0; a<n_actions; a++) {
            p[s][a] = 1.0 / (real) n_actions;
            //printf ("R(%d,%d) = %f\n", s, a, mdp->getExpectedReward(s, a));
        }
    }

    AverageValueIteration value_iteration(mdp);

    std::vector<real> Q(n_actions);
    
    value_iteration.ComputeStateValues(0.01, 1000);
    if (value_iteration.max_iter_reached) {
        for (int s=0; s<mdp->GetNStates(); s++) {
            printf ("%1.f ", value_iteration.getValue(s));
        }
        printf(" [V2] \n");
        value_iteration.ComputeStateValues(0.01, 1000);
    }
    for (int s=0; s<mdp->GetNStates(); s++) {
        printf ("%1.f ", value_iteration.getValue(s));
    }
    printf("\n");

    if (value_iteration.max_iter_reached) {
        FILE* f = fopen("test.dot", "w");
        if (f) {
            fprintf (f, "digraph MDP {\n");
            fprintf (f, "ranksep=2; rankdir=LR; \n");
            mdp->dotModel(f);
            fprintf (f, "}\n");
            fclose(f);
        }
    }

#if 0
    // GRIDOWRLD CODE
    for (int iter=0; iter < 100; iter++) {
        printf ("ITER: %d\n", iter);
        for (uint y=0; y<height; y++) {
            for (uint x=0; x<width; x++) {
                int s= x + y*width;
                for (int a=0; a<n_actions; a++) {
                    Q[a] = policy_evaluation.getValue(s,a);
                }
                //int a_max = ArgMax(Q);
                //real Q_max = Q[a_max];
                real V = policy_evaluation.getValue(s);
                Gridworld::MapElement element = grid_world.whatIs(x, y);
                if (x) printf ("&");
                if (element == Gridworld::WALL) {
                    printf ("\\bsqr ");
                } else {
                    printf ("%+.1f ", V);
                }
            }
            printf ("\\\\\\hline\n");
        }


        for (uint y=0; y<height; y++) {
            for (uint x=0; x<width; x++) {
                int s = grid_world.getState(x, y);
                Gridworld::MapElement element = grid_world.whatIs(x, y);
                if (x) printf ("&");

                if (element == Gridworld::WALL) {
                    //printf ("#");
                    printf ("\\msqr ");
                } else {
                    for (int a=0; a<n_actions; a++) {
                        Q[a] = policy_evaluation.getValue(s,a);
                    }
                    int a_max = ArgMax(Q);
                    switch (a_max) {
                    case Gridworld::NORTH: printf ("$\\uparrow$ "); break;
                    case Gridworld::SOUTH: printf ("$\\downarrow$ "); break;
                    case Gridworld::EAST:  printf ("$\\rightarrow$ "); break;
                    case Gridworld::WEST:  printf ("$\\leftarrow$ "); break;
                    }
                }
            }
            printf ("\\\\\\hline\n");
            //printf ("\n");
        }
        policy_evaluation.ComputeStateValues(0.001, 1);
        //policy_evaluation.ComputeStateActionValues(0.001, 1);
    }
    grid_world.Show();
#endif
    Stats stats = {value_iteration.Delta, value_iteration.max_iter_reached};
    return stats;
}


#endif
