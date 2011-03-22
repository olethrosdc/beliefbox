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

#include "PolicyIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

PolicyIteration::PolicyIteration(PolicyEvaluation* evaluation_,
                                 const DiscreteMDP* mdp_, 
                                 real gamma_,
                                 real baseline_) 
    : evaluation(evaluation_), mdp(mdp_), gamma(gamma_), baseline(baseline_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    
    n_actions = mdp->GetNActions();
    n_states = mdp->GetNStates();
    
    policy = new FixedDiscretePolicy(n_states, n_actions);
    evaluation->policy = policy;
    _evaluation = NULL;
    
    a_max.resize(n_states);
    
    Reset();
}

PolicyIteration::PolicyIteration(const DiscreteMDP* mdp_, 
                                 real gamma_,
                                 real baseline_) 
    : mdp(mdp_), gamma(gamma_), baseline(baseline_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    
    n_actions = mdp->GetNActions();
    n_states = mdp->GetNStates();
    
    policy = new FixedDiscretePolicy(n_states, n_actions);
    _evaluation = new PolicyEvaluation(policy, mdp, gamma, baseline);
    evaluation = _evaluation;
    
    a_max.resize(n_states);
    
    Reset();
}

void PolicyIteration::Reset()
{
    policy->Reset();
    for (int s=0; s<n_states; s++) {
        a_max[s] = ArgMax(policy->getActionProbabilitiesPtr(s));
    }
}

PolicyIteration::~PolicyIteration()
{
    delete _evaluation;
    delete policy;
}

/** ComputeStateValues
   
    threshold - exit policy estimation when difference in Q is smaller than the threshold
    max_iter - exit policy estimation when the number of iterations reaches max_iter
    The policy iteration itself stops whenever the policy has been the same for a bit.

*/
void PolicyIteration::ComputeStateValues(real threshold, int max_iter)
{
    bool policy_stable = true;
    int iter = 0;
    do {
        policy_stable = true;
        Delta = 0.0;
        // evaluate policy
        evaluation->ComputeStateValues(threshold, max_iter);
        // improve policy
        for (int s=0; s<n_states; s++) {
            real max_Qa = evaluation->getValue(s, 0);
            int argmax_Qa = 0;
            for (int a=1; a<n_actions; a++) {
                real Qa = evaluation->getValue(s, a);
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
            if (a_max[s] != argmax_Qa) {
                policy_stable = false;
                a_max[s] = argmax_Qa;
            }
        }
        Delta = evaluation->Delta;
        baseline = evaluation->baseline;
        if (max_iter >= 0) {
            iter++;
        }
    } while(policy_stable == false && iter < max_iter);

    printf ("iter left = %d\n", max_iter - iter);
}



