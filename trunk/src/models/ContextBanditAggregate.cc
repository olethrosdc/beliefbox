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
ContextBanditAggregate::ContextBanditAggregate (int n_aggregated_states, int n_states, int n_actions, real tau, real mu_0, real tau_0) :
    ContextBanditGaussian(n_states, n_actions, tau, mu_0, tau_0)

{
    mdp_dbg("Creating ContextBanditAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);

    BuildRandomAggregate();
}





/// This constructor makes a random set 
ContextBanditAggregate::ContextBanditAggregate (bool blah, int n_factors, int n_values, int n_aggregated_states, int n_states, int n_actions, real tau, real mu_0, real tau_0) :
    ContextBanditGaussian(n_states, n_actions, tau, mu_0, tau_0)

{
    mdp_dbg("Creating ContextBanditAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);

    BuildFactoredAggregate(n_factors, n_values);
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
        int s = (int) floor(urandom() * ((real) n_states));
        while (zeros && X[s].size()!=0) {
            s = (s+1) % n_states;
        }
        SMART_ASSERT(s>=0 && s<n_states);
            
        X[s].add(i);
        state_map[i] = s;
    }

}

void ContextBanditAggregate::BuildFactoredAggregate(int n_factors, int n_values)
{
    for (int i=0; i<n_aggregated_states; i++) {
        int s = i  % n_states;
        X[s].add(i);
        //printf ("Putting state %d in aggregate %d\n", i, s);
        state_map[i] = s;
    }
    int zeros = 0;
    for (int x=0; x<n_states; x++) {
        if (X[x].size()==0) {
            zeros++;
        }
    }
    if (zeros) {
        Serror ("There should not be any zeros! Found %d\n", zeros);
        exit(-1);
    }
}

ContextBanditAggregate::~ContextBanditAggregate()
{
#if 0
    printf ("Aggregate state values:\n");
    for (int s=0; s<n_aggregated_states; s++) {
        for (int a=0; a<n_actions; ++a) {
            printf ("Q(%d, %d) = %f\n", 
                    s,
                    a,
                    getExpectedReward(s,a));
        }
    }
#endif
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
    return 1.0 / (real) n_states;
}

Vector ContextBanditAggregate::getTransitionProbabilities (int s, int a) const
{
     Vector P(n_states);
    for (int i=0; i<n_states; i++) {
        P[i] = 1.0 / (real) n_states;
    }
    return P;
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
