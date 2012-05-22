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

#include "AveragePolicyEvaluation.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

AveragePolicyEvaluation::AveragePolicyEvaluation(FixedDiscretePolicy* policy_,
                                                 const DiscreteMDP* mdp_, 
                                                 real baseline_) 
    : PolicyEvaluation(policy_, mdp_, 1.0, baseline_)
{
}


AveragePolicyEvaluation::~AveragePolicyEvaluation()
{
}

/** ComputeStateValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/
void AveragePolicyEvaluation::ComputeStateValues(real threshold, int max_iter)
{
    std::vector<real> dV(n_states);
    std::vector<real> pV(n_states);
    for (int s=0; s<n_states; s++) {
        dV[s] = 0.0;
    }
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
            pV[s] =0.0;
            //printf ("V[%d] = %f, V' = ", s, V[s]);
            for (int a=0; a<n_actions; a++) {
                real p_sa = policy->getActionProbability(s, a);
                real V_sa = getValue(s, a);
                //printf ("+ %f*%f ", p_sa, V_sa);
                pV[s] += p_sa * V_sa;
            }
            //printf (" =  %f\n", pV[s]);
            dV[s] = V[s] - pV[s];
            V[s] = pV[s];
        }
	
        Delta = Max(n_states, &dV[0]) - Min(n_states, &dV[0]);
			
        if (max_iter > 0) {
            max_iter--;
        }
    } while((Delta >= threshold)  && max_iter);
    
}


