// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ContextBanditCollection.h"
#include "MultinomialDistribution.h"

ContextBanditCollection::ContextBanditCollection(int n_aggregates, int n_states, int n_actions, real tau, real mu_0, real tau_0) : MDPModel(n_states, n_actions)
{
    mdp_dbg("Creating ContextBanditCollection of %d aggregates, from %d states and %d actions\n",  n_aggregates, n_states, n_actions);
    P.resize(n_aggregates);
    A.resize(n_aggregates);
    for (int i=0; i<n_aggregates; i++) {
        P[i] = 1.0 / (real) n_aggregates;
        if (i==0) {
            A[i] = new ContextBanditGaussian(n_states, n_actions, tau, mu_0, tau_0);
        } else {
            int n = (int) floor(log(n_states) / log(2));
            // make sure 1 < n_sets < n_states
            int n_sets = 4;//1 + n_states >> (((i-1) % n) + 1);
            mdp_dbg ("Collection size: %d / %d\n", n_sets, n_states);
            A[i] = new ContextBanditAggregate(false, 3, 2, n_states, n_sets, n_actions, tau, mu_0, tau_0);
        }
    }
}


ContextBanditCollection::~ContextBanditCollection()
{

    for (uint i=0; i<A.size(); ++i) {
        printf (" %f", P[i]);
    }
    printf("# Final model probabilities\n");

    for (uint i=0; i<A.size(); i++) {
        delete A[i];
    }
}

void ContextBanditCollection::AddTransition(int s, int a, real r, int s2)
{
    // update top-level model
    real Z = 0.0;
    //printf (" L = ");
    for (uint i=0; i<A.size(); ++i) {
        //real Pi =  P[i];
        //real Psi = A[i]->getTransitionProbability(s, a, s2);
        real Psi = A[i]->getRewardDensity(s, a, r);
        //printf (" %f ", Psi);
        P[i] *= Psi;
        Z += P[i]; //post[i];
    }
    //printf("\n");
    // normalise
    if (Z > 0.0) {
        real invZ = 1.0 / Z;
        for (uint i=0; i<A.size(); ++i) {
            P[i] *= invZ;
        }
    } else {
        Swarning("normalisation failed\n");
        for (uint i=0; i<A.size(); ++i) {
            P[i] = 1.0 / (real) A.size();
        }
    }

    // update models
    //printf ("P =");
    for (uint i=0; i<A.size(); ++i) {
        //printf (" %f ", P[i]);
        A[i]->AddTransition(s, a, r, s2);
    }
    //    printf("\n");
}
real ContextBanditCollection::GenerateReward (int s, int a) const
{
    int i = DiscreteDistribution::generate(P);
    return A[i]->GenerateReward(s, a);
}
int ContextBanditCollection::GenerateTransition (int s, int a) const
{
    int i = DiscreteDistribution::generate(P);
    return A[i]->GenerateTransition(s, a);
}
real ContextBanditCollection::getTransitionProbability (int s, int a, int s2) const
{
    return 1.0 / (real) n_states;
}

real ContextBanditCollection::getRewardDensity (int s, int a, real r) const
{
    real p = 0.0;
    for (uint i=0; i<A.size(); ++i) {
        p += P[i]*A[i]->getRewardDensity(s, a, r);
    }
    return p; 
}

real ContextBanditCollection::getExpectedReward (int s, int a) const
{
#if 0
    real E = 0.0;
    for (uint i=0; i<A.size(); ++i) {
        //printf ("%f ", P[i]);
        E += P[i]*A[i]->getExpectedReward(s, a);
    }
    //printf ("= %f\n", E);
    return E; 
#else
    int i = ArgMax(P);
    return A[i]->getExpectedReward(s, a);
#endif

}
void ContextBanditCollection::Reset()
{
    //printf ("P =");
    for (uint i=0; i<A.size(); ++i) {
        //printf (" %f", P[i]);
        A[i]->Reset();
    }  
    //printf("\n");
}
