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

#include "ContextBanditGaussian.h"
#include "Random.h"

ContextBanditGaussian::ContextBanditGaussian (int n_states, int n_actions, real tau, real mu_0, real tau_0) :
    MDPModel(n_states, n_actions)
{
    mdp_dbg("Creating ContextBanditGaussian with %d states and %d actions\n",  n_states, n_actions);
    N = n_states * n_actions;
    ER.resize(N);
    for (int i=0; i<N; ++i) {
        ER[i].tau = tau;
        ER[i].mu_0 = mu_0;
        ER[i].tau_0 = tau_0;
        ER[i].Reset();
    }
}

ContextBanditGaussian::~ContextBanditGaussian()
{
    printf ("CBG state values:\n");
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; ++a) {
            printf ("Q(%d, %d) = %f\n", 
                    s,
                    a,
                    getExpectedReward(s,a));
        }
    }
    
    //printf ("COUNTS MODEL\n");
    //ShowModel();
}

void ContextBanditGaussian::AddTransition(int s, int a, real r, int s2)
{
    int ID = getID (s, a);
    //printf ("(%d, %d) [%.2f] -> %d\n", s, a, r, s2);
    ER[ID].calculatePosterior(r);
}

//void ContextBanditGaussian::SetNextReward(int s, int a, real r)
//{
//    ER[getID (s, a)].mean = r;
//}

real ContextBanditGaussian::GenerateReward (int s, int a) const
{
    return ER[getID (s, a)].generate();
}

int ContextBanditGaussian::GenerateTransition (int s, int a) const
{
    return rand()%n_states;
}

real ContextBanditGaussian::getTransitionProbability (int s, int a, int s2) const
{
    return 1.0 / (real) n_states;
}

real ContextBanditGaussian::getRewardDensity(int s, int a, real r) const
{
    return ER[getID (s,a)].pdf(r);
}

Vector ContextBanditGaussian::getTransitionProbabilities (int s, int a) const
{
    Vector P(n_states);
    for (int i=0; i<n_states; i++) {
        P[i] = 1.0 / (real) n_states;
    }
    return P;
}

real ContextBanditGaussian::getExpectedReward (int s, int a) const
{
    return ER[getID (s,a)].getMean();
}

void ContextBanditGaussian::Reset()
{
}


void ContextBanditGaussian::ShowModel() const
{
   for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "R(" << a << "," << i 
                      << ") = " << getExpectedReward(i, a)
                      << " [" << ER[getID(i,a)].n << "]"
                      << std::endl; 
        }
   }
}
