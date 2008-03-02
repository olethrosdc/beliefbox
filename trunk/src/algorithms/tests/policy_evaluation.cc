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
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"

int main (void)
{
    int period = 30;
    int max_items = 10;
    real gamma = 0.9;
    real demand = 0.1;
    InventoryManagement inventory_management (period, max_items, demand);
    
    const DiscreteMDP* mdp = inventory_management.getMDP();
    
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

    policy_evaluation.ComputeStateValues(0.01, 1000);
    for (int s=0; s<n_states; s++) {
        std::cout << s << " : " << getValue(s) << std::endl;
    }
    
}


#endif
