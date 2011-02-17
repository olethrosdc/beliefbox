// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: PolicyEstimation.c,v 1.5 2006/11/08 17:20:17 cdimitrakakis Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "PolicyEvaluation.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

PolicyEvaluation::PolicyEvaluation(DiscretePolicy* policy_,
                                   const DiscreteMDP* mdp_, 
                                   real gamma_,
                                   real baseline_) 
    : policy(policy_), mdp(mdp_), gamma(gamma_), baseline(baseline_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);

    n_actions = mdp->GetNActions();
    n_states = mdp->GetNStates();
    Reset();
}

void PolicyEvaluation::Reset()
{
    if (policy) {
        policy->Reset();
    }
    V.resize(n_states);
    for (int s=0; s<n_states; s++) {
        V[s] = 0.0;
	
    }
}

PolicyEvaluation::~PolicyEvaluation()
{
}

/** ComputeStateValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/
void PolicyEvaluation::ComputeStateValues(real threshold, int max_iter)
{
    
    int n_iter = 0;
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
            real pV =0.0;
            //printf ("S: %d ", s);
            for (int a=0; a<n_actions; a++) {
                real p_sa = policy->getActionProbability(s, a);
                real V_sa = getValue(s, a);
                //printf ("+ %f*%f \n", p_sa, V_sa);
                pV += p_sa * V_sa;
            }
            //printf (" =  %f\n", pV[s]);
            Delta += fabs(V[s] - pV);
            V[s] = pV;
        }
        
        if (max_iter > 0) {
            max_iter--;
        }
        n_iter++;
    } while((Delta >= threshold)  && max_iter != 0);
    //printf ("Exiting at delta = %f, after %d iter\n", Delta, n_iter);
}

real PolicyEvaluation::getValue (int state, int action) const
{
    real S = 0.0;
    //for (int s2=0; s2<n_states; s2++) {
    DiscreteStateSet next = mdp->getNextStates(state, action);
    for (DiscreteStateSet::iterator i=next.begin();
         i!=next.end();
         ++i) {
        int s2 = *i;
        real P = mdp->getTransitionProbability(state, action, s2);
        real R = mdp->getExpectedReward(state, action) + gamma*V[s2]  - baseline;
        S += P*R;
    }
    return S;
}



