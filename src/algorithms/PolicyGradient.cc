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

#include "PolicyGradient.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

PolicyGradient::PolicyGradient(PolicyEvaluation* evaluation_,
							   const DiscreteMDP* mdp_,
							   const Vector& starting_,
							   real gamma_)
    : evaluation(evaluation_), mdp(mdp_), starting(starting_), gamma(gamma_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    
    n_actions = mdp->getNActions();
    n_states = mdp->getNStates();
    
    policy = new FixedDiscretePolicy(n_states, n_actions);
    evaluation->policy = policy;
    
    a_max.resize(n_states);
    
    Reset();
}


void PolicyGradient::Reset()
{
    policy->Reset();
    for (int s=0; s<n_states; s++) {
        a_max[s] = ArgMax(policy->getActionProbabilitiesPtr(s));
    }
}

PolicyGradient::~PolicyGradient()
{
    delete policy;
}

/** ComputeStateValues
   
    threshold - exit policy gradient is smaller than a threshold
    max_iter - exit when the number of iterations reaches max_iter
    The policy iteration itself stops whenever the policy has been the same for a bit.

*/
void PolicyGradient::ComputeStateValues(real threshold, int max_iter)
{
    bool policy_stable = true;
    int iter = 0;	
	//Matrix G (n_states, n_actions);
    do {
        policy_stable = false;
        Delta = 0.0;
        // evaluate policy
        evaluation->ComputeStateValuesFeatureExpectation(threshold, max_iter);
		Matrix& F = evaluation->FeatureMatrix;

        if (max_iter >= 0) {
            iter++;
        }
    } while(policy_stable == false && iter < max_iter);

    printf ("iter left = %d\n", max_iter - iter);
}



