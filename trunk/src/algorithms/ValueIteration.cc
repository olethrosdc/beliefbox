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
            dQ(s, a) = 1.0;
            pQ(s, a) = 0.0;
        }
    }
}

ValueIteration::~ValueIteration()
{
}

/** Compute state values using value iteration.

	The process ends either when the error is below the given threshold,
	or when the given number of max_iter iterations is reached. Setting
	max_iter to -1 means there is no limit to the number of iterations.
*/
void ValueIteration::ComputeStateValuesStandard(real threshold, int max_iter)
{
    int n_iter = 0;
    do {
        Delta = 0.0;
        pV = V;
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
    //printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", Delta, n_iter);
}

/** Compute state values using value iteration with action elimination.

	The process ends either when the error is below the given threshold,
	or when the given number of max_iter iterations is reached. Setting
	max_iter to -1 means there is no limit to the number of iterations.

    The procedure is based on Corollary 6.7.5. in Puterman's "Markov
    Decision Proceses".

    For each iterate, we calculate the gain for state \f$s\f$:
    \f[
    B(s) = \max_a Q(s,a) - V(s) = V'(s) - V(s).
    \f]
    If, for some state-action pair \f$s,a\f$, it holds that:
    \f[
    \frac{\gamma}{1 - \gamma} \mathrm{span}(B) 
    <
    V'(s) - Q(s,a)
    \f]
    then action \f$a\f$ is sub-optimal for state \f$s\f$.
    
    \bug This does not seem to match the standard value iteration.
*/
void ValueIteration::ComputeStateValuesElimination(real threshold, int max_iter)
{
    int n_iter = 0;
    do {
        Delta = 0.0;
        pV = V;
        pQ = Q;
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                if (dQ(s,a) < 0) continue;
                real Q_sa = 0.0;
                const DiscreteStateSet& next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) - baseline;
                    Q_sa += P * (R + gamma * pV(s2));
                }
                Q(s, a) = Q_sa;
            }
            V(s) = Max(Q.getRow(s));
            dV(s) = V(s) - pV(s);
            Delta += fabs(dV(s));
        }
        
        real scale = Span(dV) * gamma / (1.0 - gamma);
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                if (dQ(s,a) < 0) continue;
                dQ(s,a) = scale  + Q(s,a) - V(s);
                //if (dQ(s,a) < 0) {
                //printf ("State %d: eliminated action %d\n", s, a);
                //}
            }
        }
        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while(Delta > threshold && max_iter != 0);
    //printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", Delta, n_iter);
}


/** Compute state values using asynchronous value iteration.

	The process ends either when the error is below the given threshold,
	or when the given number of max_iter iterations is reached. Setting
	max_iter to -1 means there is no limit to the number of iterations.

    This version updates the current values immediately
*/
void ValueIteration::ComputeStateValuesAsynchronous(real threshold, int max_iter)
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
    //printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", Delta, n_iter);
}




/** Compute state-action values using value iteration.

	The process ends either when the error is below the given threshold,
	or when the given number of max_iter iterations is reached. Setting
	max_iter to -1 means there is no limit to the number of iterations.
*/
void ValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
    ComputeStateValuesElimination(threshold, max_iter);
}

/// Create the greedy policy with respect to the calculated value function.
FixedDiscretePolicy* ValueIteration::getPolicy()
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
