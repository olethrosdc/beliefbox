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
    mdp_dbg("Creating DiscreteMDPCounts with %d states and %d actions\n",  n_states, n_actions);
    N = n_states * n_actions;
    P.resize(N);
    ER.resize(N);
    for (int i=0; i<N; ++i) {
        P[i].resize(n_states, (real) init_transition_count);
        ER[i].reset(init_reward, init_reward_count);
    }
}

DiscreteMDPCounts::~DiscreteMDPCounts()
{
    //printf ("COUNTS MODEL\n");
    //ShowModel();
}

void DiscreteMDPCounts::AddTransition(int s, int a, real r, int s2)
{
    int ID = getID (s, a);
    //printf ("(%d, %d) [%.2f] -> %d\n", s, a, r, s2);
    P[ID].Observation(s2);
    ER[ID].Observation(r);
}

//void DiscreteMDPCounts::SetNextReward(int s, int a, real r)
//{
//    ER[getID (s, a)].mean = r;
//}

real DiscreteMDPCounts::GenerateReward (int s, int a) const
{
    return ER[getID (s, a)].Generate();
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
    return ER[getID (s,a)].GetMean();
}

void DiscreteMDPCounts::Reset()
{
}


void DiscreteMDPCounts::ShowModel() const
{
    for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << a << "," << i << ":";
            for (int j=0; j<n_states; j++) {
                real p = getTransitionProbability(i, a, j);
                if (p<0.01) p =0.0f;
                std::cout << p << " ";
            }
            std::cout << " ["
                      << P[getID(i,a)].GetParameters().Sum()
                      << "]\n";
        }
    }

   for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "R(" << a << "," << i 
                      << ") = " << getExpectedReward(i, a)
                      << " [" << ER[getID(i,a)].n_samples << "]"
                      << std::endl; 
        }
   }
}
DiscreteMDP* DiscreteMDPCounts::generate() 
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions, NULL, NULL);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Vector C =  P[getID (s,a)].generate();
            real expected_reward = getExpectedReward(s,a);
            mdp->addFixedReward(s, a, expected_reward);
            for (int s2=0; s2<n_states; s2++) {
                mdp->setTransitionProbability(s, a, s2, C[s2]);
            }
        }
    }
    
    return mdp;
}


DiscreteMDP* DiscreteMDPCounts::getMeanMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions, NULL, NULL);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Vector C =  P[getID (s,a)].GetMean();
            real expected_reward = getExpectedReward(s,a);
            mdp->addFixedReward(s, a, expected_reward);
            for (int s2=0; s2<n_states; s2++) {
                mdp->setTransitionProbability(s, a, s2, C[s2]);
            }
        }
    }
    
    return mdp;
}
