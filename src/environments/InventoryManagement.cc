// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "InventoryManagement.h"
#include "SpecialFunctions.h"



//  use ignbin() for generation.

InventoryManagement::InventoryManagement(int period_,
                                         int max_items_,
                                         real demand_,
                                         real margin_)
    : period(period_), max_items(max_items_), demand(demand_), margin(margin_)
{
    n_actions = max_items + 1;
    n_states = max_items + 1;
    
    assert (max_items > 1);
    assert (margin >= 1);
    assert (demand >= 0 && demand <= 1);
    assert (period > 0);
    local_mdp = getMDP();
    local_mdp->Check();
    Reset();
}

DiscreteMDP* InventoryManagement::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions, NULL);
#if 0
    for (uint s=0; s<n_states; s++) {
        for (uint a=0; a<n_actions; a++) {
            int order = a;
            if (s+a > (uint) max_items) {
                order = max_items - s;
            }
            real expected_reward = - order; // add the order costs first

            // Calculate the transition probabilities
            int current_stock = s + order;
            real P_empty = 1;
            for (int s2 = 0; s2<n_states; ++s2) {
                mdp->setTransitionProbability(s, a, s2, 0.0);
            }
            for (int d=0; d<current_stock; ++d) {
                int s2 = current_stock - d;
                real P =  binomial(period, d)*pow(demand, d)*pow(1-demand, period - d);
                P_empty -= P;
                expected_reward += margin * P*d;
                mdp->setTransitionProbability(s, a, s2, P);
            }
            // Special case for empty stock due to large demand.
            expected_reward += margin * P_empty * current_stock;
            mdp->setTransitionProbability(s, a, 0, P_empty);

            // Add the expected reward.
            mdp->setFixedReward(s, a, expected_reward);
        }
    }
#else
    for (uint s=0; s<n_states; s++) {
        for (uint a=0; a<n_actions; a++) {
            mdp->setFixedReward(s, a, (real) s + 0.1 * (real) a);
            Vector p(n_states);
            for (uint s2 = 0; s2<n_states; ++s2) {
                p(s2) = urandom();
            }
            p /= p.Sum();
            for (uint s2 = 0; s2<n_states; ++s2) {
                mdp->setTransitionProbability(s, a, s2, p(s2));
            }
        }
        
    }
#endif
    mdp->Check();
    mdp->ShowModel();
    return mdp;
}

InventoryManagement::~InventoryManagement()
{
    delete local_mdp;
}
