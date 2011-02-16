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

#include "ValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

ValueIteration::ValueIteration(const DiscreteMDP* mdp, real gamma, real baseline)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    this->mdp = mdp;
    this->gamma = gamma;
    this->baseline = baseline;
    n_actions = mdp->GetNActions();
    n_states = mdp->GetNStates();
    Reset();
}

void ValueIteration::Reset()
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
            dQ(s, a) = 0.0;
            pQ(s, a) = 0.0;
        }
    }
}

ValueIteration::~ValueIteration()
{
}

void ValueIteration::ComputeStateValues(real threshold, int max_iter)
{
        //Vector pV(V.size());
        //Vector dV(V.size());
    
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
                //real v = V(s);
            real Q_a_max = -RAND_MAX;
            int a_max = 0;
            for (int a=0; a<n_actions; a++) {
                real Q_sa = 0.0;
                const DiscreteStateSet& next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) + gamma * V(s2) - baseline;
                    Q_sa += P * R;
                }
                Q(s, a) = Q_sa;
                if (a==0 || Q_a_max < Q_sa) {
                    a_max = a;
                    Q_a_max = Q_sa;
                }
            }
            V(s) = Q_a_max;
            dV(s) = fabs(pV(s) - V(s));
            pV(s) = V(s);
        }
        Delta = dV.Norm(1.0);
        if (max_iter > 0) {
            max_iter--;
        }
    } while(Delta >= threshold && max_iter != 0);
	
}




/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void ValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
    int N = n_states * n_actions;

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            dQ(s, a) = 0.0;
        }
    }
    int n_iter = 0;
    do {
        Delta = 0.0;
        for (int s0=0; s0<n_states; s0++) {
            int s = s0;
            for (int a=0; a<n_actions; a++) {
                real sum = 0.0;

                const DiscreteStateSet& next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) - baseline;
                    real Q_a_max = Max(Q.getRow(s2));
                    sum += P*(R + gamma*Q_a_max);
                }
                Q(s, a) = sum;
                dQ(s, a) = pQ(s, a) - sum;
                Delta += fabs(dQ(s, a));
                pQ(s, a) = sum;
            }
        }
        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    //printf("Exiting at d:%f, n:%d\n", Delta, n_iter);
}


