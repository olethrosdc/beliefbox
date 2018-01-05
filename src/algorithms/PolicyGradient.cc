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

PolicyGradient::PolicyGradient(const DiscreteMDP* mdp_,
                               real gamma_,
                               real step_size_)
    : evaluation(NULL, mdp_, gamma_),
      mdp(mdp_),
      starting(mdp->getNStates()),
      gamma(gamma_),
      step_size(step_size_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    
    n_actions = mdp->getNActions();
    n_states = mdp->getNStates();
    
    policy = new FixedDiscretePolicy(n_states, n_actions);
    evaluation.policy = policy;
    
    a_max.resize(n_states);
    for (int i=0; i<n_states; ++i) {
        starting(i) = 1.0 / (real) n_states;
    }
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

*/
void PolicyGradient::ModelBasedGradient(real threshold, int max_iter)
{
    bool policy_stable = true;
    int iter = 0;	
	//Matrix G (n_states, n_actions);
    do {
        policy_stable = false;
        Delta = 0.0;
        // evaluate policy
        evaluation.ComputeStateValuesFeatureExpectation(threshold, max_iter);
        Matrix& F = evaluation.FeatureMatrix;

        // r = E[reward | s]
        Vector r(n_states);
        for (int i=0; i<n_states; i++) {
            for (int a=0; a<n_actions; a++) {
                r(i) += mdp->getExpectedReward(i, a) * policy->getActionProbability(i, a);
            }
        }
        // L = y' X
        Vector L(n_states);
        for (int i=0; i<n_states; i++) {
            for (int j=0; j<n_states; j++) {
                L(i) += starting(j) * F(i, j);
            }
        }

        // R = X r
        Vector R(n_states);
        for (int i=0; i<n_states; i++) {
            for (int j=0; j<n_states; j++) {
                R(i) += F(i, j) * r(j);
            }
        }


        Matrix D(n_states, n_actions);
        for (int a=0; a<n_actions; a++) {
            for (int i=0; i<n_states; i++) {
                for (int j=0; j<n_states; j++) {
                    D(i, a) += mdp->getTransitionProbability(i, a, j) * R(j);
                }
                D(i, a) *= L(i);
                Delta += fabs(D(i,a));
            }
        }

        for (int i=0; i<n_states; i++) {
            Vector* Theta = policy->getActionProbabilitiesPtr(i);
            for (int a=0; a<n_actions; a++) {
                (*Theta)(a) += step_size * D(i, a);
            }
        }

        if (Delta < 0.1) {
            policy_stable = true;
        }

        if (max_iter >= 0) {
            iter++;
        }
    } while(policy_stable == false && iter < max_iter);

    printf ("iter left = %d\n", max_iter - iter);
}



