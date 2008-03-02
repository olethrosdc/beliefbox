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
    policy->Reset();
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
    std::vector<real> dV(n_states);
    std::vector<real> pV(n_states);
    for (int s=0; s<n_states; s++) {
        dV[s] = 0.0;
    }
    do {
        Delta = 0.0;
        for (int s0=0; s0<n_states; s0++) {
            int s = s0;
	    pV[s] =0.0;
            for (int a=0; a<n_actions; a++) {
		pV[s] += policy->getActionProbability(s, a) * getValue(s, a);
            }
	    dV[s] = V[s] - pV[s];
        }
	for (int s=0; s<n_states; s++) {
	    V[s] = pV[s];
	}
		
	Delta = Max(n_states, &dV[0]) - Min(n_states, &dV[0]);
			
        if (max_iter > 0) {
            max_iter--;
            //	printf ("left: %d\n", max_iter);
        }
    } while((Delta >= threshold)  && max_iter);
    
    //if (!max_iter) {
    //   fprintf (stderr, "warning - delta %f >= %f\n", Delta, threshold);
    //}		
}

real PolicyEvaluation::getValue (int state, int action)
{
    real sum_over_states = 0.0;
    for (int s2=0; s2<n_states; s2++) {
	real P = mdp->getTransitionProbability(state, action, s2);
	real R = mdp->getExpectedReward(state, action) - baseline;
	sum_over_states += P*(R + gamma*V[s2]);
    }
    return sum_over_states;
}


real PolicyEvaluation::getValue (int state)
{
    return V[state];
}
