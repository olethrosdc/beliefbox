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

real DiscreteMDPAggregate::GenerateReward (int s, int a)
{
    return DiscreteMDPCounts::GenerateReward(Aggregate(s), a);
}

int DiscreteMDPAggregate::GenerateTransition (int s, int a)
{
    return DiscreteMDPCounts::GenerateTransition(Aggregate(s), a);
}

real DiscreteMDPAggregate::getTransitionProbability (int s, int a, int s2)
{
    return DiscreteMDPCounts::getTransitionProbability(Aggregate(s), a, Aggregate(s2));
}

Vector DiscreteMDPAggregate::getTransitionProbabilities (int s, int a)
{
    return DiscreteMDPCounts::getTransitionProbabilities(Aggregate(s), a);
}

real DiscreteMDPAggregate::getExpectedReward (int s, int a)
{
    return DiscreteMDPCounts::getExpectedReward(Aggregate(s), a);
}

void DiscreteMDPAggregate::Reset()
{
    DiscreteMDPAggregate::Reset();
}
