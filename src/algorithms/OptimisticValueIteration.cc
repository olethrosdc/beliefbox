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
#include "ValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include "DiscretePolicy.h"
#include "Bounds.h"
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
    V.Resize(n_states);
    V.Clear();
    Q.Resize(n_states, n_actions);
    Q.Clear();
}

OptimisticValueIteration::~OptimisticValueIteration()
{
}

/** Compute state values using value iteration in an augmented MDP.

    The augmented MDP is 

    \param delta the error probability for the confidence bound.
    \param threshold stop when the change in value Delta < threshold.
    \param max_iter stop after at most max_iter steps, unless max_iter < 0.
*/
void OptimisticValueIteration::ComputeStateValuesAugmentedMDP(real delta,
                                                              real threshold,
                                                              int max_iter)
{
    int n_aug = 2 * n_actions * n_states;
    DiscreteMDP augmented_mdp(n_states, n_aug);
    for (int s=0; s<n_states; s++) {
        int a_aug = 0;
        for (int a=0; a<n_actions; a++) {
            int N_sa = (int) ceil(mdp->getNVisits(s, a));
            real r_sa = mdp->getExpectedReward(s,a);
            real r_gap = HoeffdingBound(1, N_sa, delta);
            Vector P_sa = mdp->getTransitionProbabilities(s, a);
            real gap = WeissmanBound(n_states, N_sa, delta);
            for (int k=0; k<n_states; k++) {
                augmented_mdp.setTransitionProbabilities(s, a_aug, MultinomialDeviation(P_sa, k, -gap));
                augmented_mdp.setFixedReward(s, a_aug, r_sa + r_gap);
                a_aug++;
                augmented_mdp.setTransitionProbabilities(s, a_aug, MultinomialDeviation(P_sa, k, +gap));
                augmented_mdp.setFixedReward(s, a_aug, r_sa + r_gap);
                a_aug++;
            }
        }
    }

    augmented_mdp.Check();

    ValueIteration vi(&augmented_mdp, gamma);
    vi.ComputeStateValues(threshold, max_iter);
    for (int s=0; s<n_states; s++) {
        int a_aug = 0;
        V(s) = -INF;
        for (int a=0; a<n_actions; a++) {
            real Qaug_max = vi.getValue(s, a_aug);
            for (int k=0; k<2*n_states; ++k, a_aug++) {
                //printf("%d %f\n", k, vi.getValue(s, a_aug));
                Qaug_max = std::max(Qaug_max, vi.getValue(s, a_aug));
            }
            Q(s,a) = Qaug_max;
            V(s) = std::max(V(s), Q(s,a));
            //printf("Q(%d, %d) = %f\n", s, a, Q(s,a));
        }
    }            
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
