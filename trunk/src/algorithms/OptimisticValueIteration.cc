// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: ValueIteration.c,v 1.5 2006/11/08 17:20:17 cdimitrakakis Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "OptimisticValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include "DiscretePolicy.h"
#include <cmath>
#include <cassert>

OptimisticValueIteration::OptimisticValueIteration(const DiscreteMDPCounts* mdp,
                                                   real gamma,
                                                   real baseline)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    this->mdp = mdp;
    this->gamma = gamma;
    this->baseline = baseline;
    n_actions = mdp->getNActions();
    n_states = mdp->getNStates();
    Reset();
}

void OptimisticValueIteration::Reset()
{
    //int N = n_states * n_actions;

    V.Resize(n_states);
    dV.Resize(n_states);
    pV.Resize(n_states);

    Q.Resize(n_states, n_actions);
    dQ.Resize(n_states, n_actions);
    pQ.Resize(n_states, n_actions);

    
    for (int s=0; s<n_states; s++) {
        V(s) = 0.0;
        dV(s) = 0.0;
        pV(s) = 0.0;
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = 0.0;
            dQ(s, a) = 1.0;
            pQ(s, a) = 0.0;
        }
    }
}

OptimisticValueIteration::~OptimisticValueIteration()
{
}

/** Compute state values using value iteration.

	The process ends either when the error is below the given threshold,
	or when the given number of max_iter iterations is reached. Setting
	max_iter to -1 means there is no limit to the number of iterations.
*/
void OptimisticValueIteration::ComputeStateValuesStandard(real error_probability, real epsilon, real threshold, int max_iter)
{
    int n_iter = 0;
    do {
        Delta = 0.0;
        pV = V;
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                real Q_sa = 0.0;
                int N_sa = (int) ceil(mdp->getNVisits(s, a));
                for (int s2=0; s2<n_states; ++s2) {
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) - baseline;
                    Q_sa += P * (R + gamma * pV(s2));
                }
                Q(s, a) = Q_sa;
            }
            V(s) = Max(Q.getRow(s));
            Delta += fabs(V(s) - pV(s));
        }
        
        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    //printf("#OptimisticValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", Delta, n_iter);
}


/// Create the greedy policy with respect to the calculated value function.
FixedDiscretePolicy* OptimisticValueIteration::getPolicy() const
{
    FixedDiscretePolicy* policy = new FixedDiscretePolicy(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        int argmax_Qa = ArgMax(Q.getRow(s));
        Vector* p = policy->getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) { 
            (*p)(a) = 0.0;
        }
        (*p)(argmax_Qa) = 1.0;
    }
    return policy;
}
