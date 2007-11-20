/* -*- Mode: C++; -*- */
/* VER: $Id: VPIPolicy.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteBanditPolicy.h"


VPIPolicy::VPIPolicy(int n_actions, ActionValueEstimate* estimator, real gamma, int n_samples)
{
    this->n_actions = n_actions;
    this->estimator = estimator;
    this->gamma = gamma;
    this->n_samples = n_samples;
}

/// Reset
void VPIPolicy::Reset() {
    estimator->Reset();
}

void VPIPolicy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}

VPIPolicy::~VPIPolicy()
{
}

int VPIPolicy::SelectAction()
{
    real B = 1.0f / (1.0f - gamma); // the exploration factor
    int j = estimator->GetMax();
    int j2 = estimator->GetSecondMax();
    real q_j = B * estimator->GetMean(j); // return of action j
    real q_j2 = B * estimator->GetMean(j2); // return of action j2
    int a = j; // the selected action (initially greedy)
    real U_max = -1;  // updated when i == 0

    //printf ("E[%d] = %f, E[%d] = %f\n" ,j, q_j, j2, q_j2);

    // Estimate the utility of all acitons
    for (int i=0; i<n_actions; i++) {
        real V = B * estimator->GetMean(i);
        real Gain = 0.0f;

        // Integrate over gain
        for (int n=0; n<n_samples; n++) {
            // Calculate max policy from current value
            real sample_max = estimator->GetMean(j);
            //real sample_max = q_j;
            real q_i = estimator->Sample(i);
            real q_i_max;
            if (q_i > sample_max) {
                q_i_max = B*q_i;
            } else {
                q_i_max = q_i + gamma * B * sample_max;
            }

            if (i==j && q_j2 > q_i_max) {
                Gain += q_j2 - q_i_max;
            } else if (i!=j && q_i_max > q_j) {
                Gain += q_i_max - q_j;
            }
        }
        Gain /= (real) n_samples;
        //real U = Gain;
        real U = V + Gain;
        if (U > U_max || i == 0) {
            U_max = U;
            a = i;
        }
    }
    return a;
}

