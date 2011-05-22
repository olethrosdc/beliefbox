// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "GoalStateIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

GoalStateIteration::GoalStateIteration(const DiscreteMDP* mdp)
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

void GoalStateIteration::Reset()
{
    //int N = n_states * n_actions;
    V.resize(n_states);
    dV.resize(n_states);
    pV.resize(n_states);
    for (int s=0; s<n_states; s++) {
        V[s] = 0.0;
        dV[s] = 0.0;
        pV[s] = 0.0;
    }
}

GoalStateIteration::~GoalStateIteration()
{
}

// calculate the negative expected distance
void GoalStateIteration::ComputeStateValues(int goal_state, real threshold, int max_iter)
{
    assert (goal_state >=0 && goal_state < n_states);
    Reset();
    V[goal_state] = 1;
    
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
            if (s == goal_state) {
                continue;
            }
            real Q_a_max = -RAND_MAX;
            int a_max = 0;
            for (int a=0; a<n_actions; a++) {
                real S = 0.0;
                DiscreteStateSet next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = V[s2] - 1;
                    S += P * R;
                }
                if (a==0 || Q_a_max < S) {
                    a_max = a;
                    Q_a_max = S;
                }
            }
            V[s] = Q_a_max;
            dV[s] = pV[s] - V[s];
            pV[s] = V[s];
        }
        Delta = Max(dV) - Min(dV);
        max_iter--;

    } while(Delta >= threshold && max_iter > 0);
	
}

/// Compute max_s min_pi expected time to reach goal_state
real GoalStateIteration::GetMaximumDistanceFromState(int goal_state, real threshold, int max_iter)
{
    ComputeStateValues(goal_state, threshold, max_iter);
    
    real min_V = 0;
    for (int s=0; s<n_states; s++) {
        if (s==goal_state) {
            continue;
        }
        real Q_a_max = -RAND_MAX;
        int a_max = 0;
        for (int a=0; a<n_actions; a++) {
            real S = 0.0;
            DiscreteStateSet next = mdp->getNextStates(s, a);
            for (DiscreteStateSet::iterator i=next.begin();
                 i!=next.end();
                 ++i) {
                int s2 = *i;
                real P = mdp->getTransitionProbability(s, a, s2);
                real R = V[s2] - 1;
                S += P * R;
            }
            if (a==0 || Q_a_max < S) {
                a_max = a;
                Q_a_max = S;
            }
        }
        if (Q_a_max < min_V) {
            min_V = Q_a_max;
        }
    }	
    //printf ("min_V = %f\n", min_V);
    return - min_V;
}

/// Compute max_s expected time to reach any state
real GoalStateIteration::GetMDPDiameter(real threshold, int max_iter)
{
    real max_V = 0;
    for (int s=0; s<n_states; s++) {
        real V = GetMaximumDistanceFromState(s, threshold, max_iter);
        if (V > max_V) {
            max_V = V;
        }
    }	
    return max_V;
}

