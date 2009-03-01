// -*- Mode: c++ -*-
// copyright (c) 2005-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDPAggregate.h"
#include "Random.h"



/// This constructor makes a random set 
DiscreteMDPAggregate::DiscreteMDPAggregate (int n_aggregated_states, int n_states, int n_actions, int init_transition_count, int init_reward_count, real init_reward) :
    DiscreteMDPCounts(n_states, n_actions, init_transition_count,  init_reward_count,  init_reward)

{
    mdp_dbg("Creating DiscreteMDPAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);
    for (int i=0; i<n_aggregated_states; i++) {
        int s = floor(urandom() * ((real) n_states));
        SMART_ASSERT(s>=0 && s<n_states);
        X[s].add(i);
        state_map[i] = s;
    }
}

DiscreteMDPAggregate::~DiscreteMDPAggregate()
{
}

void DiscreteMDPAggregate::AddTransition(int s, int a, real r, int s2)
{
    DiscreteMDPCounts::AddTransition(Aggregate(s), a, r, Aggregate(s2));
}

void DiscreteMDPAggregate::SetNextReward(int s, int a, real r)
{
    DiscreteMDPCounts::SetNextReward(Aggregate(s), a, r);
}

real DiscreteMDPAggregate::GenerateReward (int s, int a) const
{
    return DiscreteMDPCounts::GenerateReward(Aggregate(s), a);
}

int DiscreteMDPAggregate::GenerateTransition (int s, int a) const
{
    return DiscreteMDPCounts::GenerateTransition(Aggregate(s), a);
}

real DiscreteMDPAggregate::getTransitionProbability (int s, int a, int s2) const
{
    int x = Aggregate(s);
    int x2 = Aggregate(s2);
    int n = X[x2].size();
    real p = 0.0;
    if (n > 0) {
        p = DiscreteMDPCounts::getTransitionProbability(x, a, x2) / (real) n; 
    }
    mdp_dbg ("p (%d | %d, %d) = %f\n", x2, x, a, p);
    return p;
}

Vector DiscreteMDPAggregate::getTransitionProbabilities (int s, int a) const
{
    return DiscreteMDPCounts::getTransitionProbabilities(Aggregate(s), a);
}

real DiscreteMDPAggregate::getExpectedReward (int s, int a) const
{
    return DiscreteMDPCounts::getExpectedReward(Aggregate(s), a);
}

void DiscreteMDPAggregate::Reset()
{
    DiscreteMDPAggregate::Reset();
}
