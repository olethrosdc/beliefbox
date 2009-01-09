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

DiscreteMDPCounts::DiscreteMDPCounts (int n_states, int n_actions, int init_transition_count, int init_reward_count, real init_reward) :
    MDPModel(n_states, n_actions)
{
    N = n_states * n_actions;
    P.resize(N);
    for (int i=0; i<N; ++i) {
        P[i].resize(n_states, (real) init_transition_count);
    }
}

DiscreteMDPCounts::~DiscreteMDPCounts();
void DiscreteMDPCounts::AddTransition(int s, int a, real r, int s2);
void DiscreteMDPCounts::ShowModel();
real DiscreteMDPCounts::GenerateReward (int s, int a);
int DiscreteMDPCounts::GenerateTransition (int s, int a);
real DiscreteMDPCounts::getTransitionProbability (int s, int a, int s2);
real DiscreteMDPCounts::getExpectedReward (int s, int a);
int DiscreteMDPCounts::getNStates () { return N;}
void DiscreteMDPCounts::Reset();
