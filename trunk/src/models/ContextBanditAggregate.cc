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

#include "ContextBanditAggregate.h"
#include "Random.h"
#include "Gridworld.h"


/// This constructor makes a random set 
ContextBanditAggregate::ContextBanditAggregate (int n_aggregated_states, int n_states, int n_actions, int init_transition_count, int init_reward_count, real init_reward) :
    ContextBanditGaussian(n_states, n_actions, init_transition_count,  init_reward_count,  init_reward)

{
    mdp_dbg("Creating ContextBanditAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);

    BuildRandomAggregate();
}


void ContextBanditAggregate::BuildRandomAggregate()
{
    for (int i=0; i<n_aggregated_states; i++) {
        int zeros = 0;
        for (int x=0; x<n_states; x++) {
            if (X[x].size()==0) {
                zeros++;
            }
        }
        int s = floor(urandom() * ((real) n_states));
        while (zeros && X[s].size()!=0) {
            s = (s+1) % n_states;
        }
        SMART_ASSERT(s>=0 && s<n_states);
            
        X[s].add(i);
        state_map[i] = s;
    }

}


ContextBanditAggregate::~ContextBanditAggregate()
{
}

void ContextBanditAggregate::AddTransition(int s, int a, real r, int s2)
{
    ContextBanditGaussian::AddTransition(Aggregate(s), a, r, Aggregate(s2));
}

//void ContextBanditAggregate::SetNextReward(int s, int a, real r)
//{
//    ContextBanditGaussian::SetNextReward(Aggregate(s), a, r);
//}

real ContextBanditAggregate::GenerateReward (int s, int a) const
{
    return ContextBanditGaussian::GenerateReward(Aggregate(s), a);
}

int ContextBanditAggregate::GenerateTransition (int s, int a) const
{
    return ContextBanditGaussian::GenerateTransition(Aggregate(s), a);
}

real ContextBanditAggregate::getTransitionProbability (int s, int a, int s2) const
{
    int x = Aggregate(s);
    int x2 = Aggregate(s2);
    int n = X[x2].size();
    real p = 0.0;
    if (n > 0) {
        p = ContextBanditGaussian::getTransitionProbability(x, a, x2) / (real) n; 
    }
    mdp_dbg ("\n P(s'=%d | s=%d, a=%d) = (1/%d) P(x'=%d | x=%d, a=%d) = %f\n",
             s2, s, a,
             n, x2, x, a,
             p);
    return p;
}

Vector ContextBanditAggregate::getTransitionProbabilities (int s, int a) const
{
    return ContextBanditGaussian::getTransitionProbabilities(Aggregate(s), a);
}

real ContextBanditAggregate::getRewardDensity (int s, int a, real r) const
{
    return ContextBanditGaussian::getRewardDensity(Aggregate(s), a, r);
}

real ContextBanditAggregate::getExpectedReward (int s, int a) const
{
    return ContextBanditGaussian::getExpectedReward(Aggregate(s), a);
}

void ContextBanditAggregate::Reset()
{
    ContextBanditGaussian::Reset();
}
