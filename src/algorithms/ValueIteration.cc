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
    n_actions = mdp->getNActions();
    n_states = mdp->getNStates();
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
    int n_iter = 0;
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                real Q_sa = 0.0;
                const DiscreteStateSet& next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) - baseline;
                    Q_sa += P * (R + gamma * V(s2));
                }
                Q(s, a) = Q_sa;
            }
            V(s) = Max(Q.getRow(s));
            Delta += fabs(V(s) - pV(s));
            pV(s) = V(s);
        }

        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    printf("Exiting at d:%f, n:%d\n", Delta, n_iter);

}




/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void ValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
    int n_iter = 0;
    do {
        Delta = 0.0;
        //for (int s0=0; s0<n_states; s0++) {
        //int s = s0;
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                real sum = 0.0;

                const DiscreteStateSet& next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) - baseline;
                    sum += P*(R + gamma*V(s2));
                    //real Q_a_max = Max(Q.getRow(s2));
                    //sum += P*(R + gamma*Q_a_max);
                }
                Q(s, a) = sum;
                Delta += fabs(pQ(s, a) - sum);
                pQ(s, a) = sum;
            }
            V(s) = Max(Q.getRow(s));
        }
        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    //printf("Exiting at d:%f, n:%d\n", Delta, n_iter);
}


FixedDiscretePolicy* ValueIteration::getPolicy()
{
    FixedDiscretePolicy* policy = new FixedDiscretePolicy(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        real max_Qa = getValue(s, 0);
        int argmax_Qa = 0;
        for (int a=1; a<n_actions; a++) {
            real Qa = getValue(s, a);
            if (Qa > max_Qa) {
                max_Qa = Qa;
                argmax_Qa = a;
            }
        }
        Vector* p = policy->getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) { 
            (*p)[a] = 0.0;
        }
        (*p)[argmax_Qa] = 1.0;
    }
    return policy;
}
