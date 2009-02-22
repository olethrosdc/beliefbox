// -*- Mode: c++ -*-
// copyright (c) 2005-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDPCounts.h"
#include "Random.h"

DiscreteMDPCounts::DiscreteMDPCounts (int n_states, int n_actions, int init_transition_count, int init_reward_count, real init_reward) :
    MDPModel(n_states, n_actions)
{
    N = n_states * n_actions;
    P.resize(N);
    ER.resize(N);
    for (int i=0; i<N; ++i) {
        P[i].resize(n_states, (real) init_transition_count);
        ER[i] =0.0;
    }
}

DiscreteMDPCounts::~DiscreteMDPCounts()
{
}

void DiscreteMDPCounts::AddTransition(int s, int a, real r, int s2)
{
    int ID = getID (s, a);
    P[ID].Observation(s2);
    ER[ID] = r;
}

void DiscreteMDPCounts::SetNextReward(int s, int a, real r)
{
    ER[getID (s, a)] = r;
}

real DiscreteMDPCounts::GenerateReward (int s, int a) const
{
    return ER[getID (s, a)];
}

int DiscreteMDPCounts::GenerateTransition (int s, int a) const
{
    Vector p = P[getID (s,a)].GetMean();
    real d=urandom();
    real sum = 0.0;
    int n_outcomes = p.Size();
    for (int i=0; i<n_outcomes; i++) {
        sum += p[i];
        if (d < sum) {
            return i;
        }
    }
    return rand()%n_outcomes;
}

real DiscreteMDPCounts::getTransitionProbability (int s, int a, int s2) const
{
    Vector p = P[getID (s,a)].GetMean();
    return p[s2];
}

Vector DiscreteMDPCounts::getTransitionProbabilities (int s, int a) const
{
    return P[getID (s,a)].GetMean();
}

real DiscreteMDPCounts::getExpectedReward (int s, int a) const
{
    return ER[getID (s,a)];
}

void DiscreteMDPCounts::Reset()
{
}


