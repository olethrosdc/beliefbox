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
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"

int main (void)
{
    int period = 30;
    int max_items = 10;
    real gamma = 1.0;
    real demand = 0.1;
    real random = 0.0;
    real pit = -100.0;
    uint width=8;
    uint height=8;
    InventoryManagement inventory_management (period, max_items, demand);

    Gridworld grid_world("maze1", width, height, 4, random, pit);
    //const DiscreteMDP* mdp = inventory_management.getMDP();
    const DiscreteMDP* mdp = grid_world.getMDP();
    
    int n_states = mdp->GetNStates();
    int n_actions = mdp->GetNActions();
    
    std::vector<Vector> p(n_states);
    for (int s=0; s<n_states; s++) {
        p[s].Resize(n_actions);
        for (int a=0; a<n_actions; a++) {
            p[s][a] = 1.0 / (real) n_actions;
        }
    }

    DiscretePolicy* policy = new FixedDiscretePolicy(p);

    PolicyEvaluation policy_evaluation(policy, mdp, gamma);
    //ValueIteration policy_evaluation(mdp, gamma);

    policy_evaluation.ComputeStateValues(0.001, 1000);
    std::vector<real> Q(n_actions);
#if 0
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q[a] = policy_evaluation.getValue(s,a);
        }
        int a_max = ArgMax(Q);
        std::cout << s << " : " << policy_evaluation.getValue(s)
                  << " | "
                  << a_max << ":" << Q[a_max]
                  << std::endl;
    }
#else
    for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
            int s= x + y*width;
            for (int a=0; a<n_actions; a++) {
                Q[a] = policy_evaluation.getValue(s,a);
            }
            int a_max = ArgMax(Q);
            real Q_max = Q[a_max];
            real V = policy_evaluation.getValue(s);
            printf ("%+.1f ", V);
        }
        printf ("\n");
    }
#endif
    
}


#endif
