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
                                         real demand_)
    : period(period_), max_items(max_items_), demand(demand_)
{
    int n_actions = max_items + 1;
    int n_states = max_items + 1;

    mdp = new DiscreteMDP(n_states, n_actions, NULL, NULL);

    
    dist.resize(period + 1);
    for (int i=0; i<=period; i++) {
        dist[i] = new BinomialDistribution(period, i, -1);
    }


    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int order = a;
            if (s+a > max_items) {
                order = max_items - s;
            }
            for (int s2=0; s2<n_states; s2++) {
                int d = s + order - s2;
                mdp->setTransitionProbability(s, a, s2, binomial(period, d)*pow(demand, d)*pow(1-demand, period - d));
                mdp->setRewardDistribution(s, a, dist[d]);
            }
        }
    }
}

InventoryManagement::~InventoryManagement()
{
    for (int i=0; i<=period; i++) {
        delete dist[i];
    }

}
