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
    int n_actions = max_items + 1;
    int n_states = max_items + 1;
    
    assert (max_items > 1);
    assert (margin >= 1);
    assert (demand >= 0 && demand <= 1);
    assert (period > 0);

    mdp = new DiscreteMDP(n_states, n_actions, NULL);

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int order = a;
            if (s+a > max_items) {
                order = max_items - s;
            }
            real expected_reward = - order; // add the order costs first

            // Calculate the transition probabilities
            int current_stock = s + order;
            real P_empty = 1;
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
            dist.push_back(new SingularDistribution(expected_reward));
            mdp->setRewardDistribution(s, a, dist.back());
        }
    }
    mdp->Check();
}

InventoryManagement::~InventoryManagement()
{
    for (int i=0; i<dist.size(); i++) {
        delete dist[i];
    }
}
