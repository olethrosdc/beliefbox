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

#include "OptimisticValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include "DiscreteMDPCounts.h"
#include <cmath>
#include <cassert>

OptimisticValueIteration::OptimisticValueIteration(DiscreteMDPCounts* mdp, real gamma, real baseline)
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

void OptimisticValueIteration::Reset()
{
    int N = n_states * n_actions;
    V.resize(n_states);
    dV.resize(n_states);
    pV.resize(n_states);
    Q.resize(n_states);
    Q_data.resize(N);
    dQ.resize(n_states);
    dQ_data.resize(N);
    pQ.resize(n_states);
    pQ_data.resize(N);
    for (int s=0; s<n_states; s++) {
        V[s] = 0.0;
        dV[s] = 0.0;
        pV[s] = 0.0;
        Q[s] = &Q_data[s*n_actions];
        dQ[s] = &dQ_data[s*n_actions];
        pQ[s] = &pQ_data[s*n_actions];
        for (int a=0; a<n_actions; a++) {
            Q[s][a] = 0.0;
            dQ[s][a] = 0.0;
            pQ[s][a] = 0.0;
        }
    }
}

OptimisticValueIteration::~OptimisticValueIteration()
{
}

void OptimisticValueIteration::ComputeStateValues(real epsilon, real threshold, int max_iter)
{
    do {
        Delta = 0.0;
        for (int s=0; s<n_states; s++) {
            //real v = V[s];
            real Q_a_max = -RAND_MAX;
            int a_max = 0;
            for (int a=0; a<n_actions; a++) {
                real S = 0.0;

                Vector V2(n_states);
                // store the value of next states
                for (int s2=0; s2<n_states; ++s2) {
                    V2[s2] = V[s2];
                }
                Vector Q = mdp->getTransitionProbabilities(s, a);
                real max_U = -RAND_MAX;
                
                for (int s2=0; s2<n_states; ++s2) {
                    Vector I(n_states);
                    I[s2] = 1.0;
                    I = I + (I - 1.0)/((real) (n_states -1));
                    for (int j=-1; j<=1; j+=2) {
                        Vector P(Q);
                        P += I*((real) j)*epsilon;
                        real U = Product(&P, &V2);
                        if (U > max_U) {
                            max_U = U;
                        }
                    }
                }
                        
                real R = mdp->getExpectedReward(s, a) + gamma*max_U - baseline;
                
                if (a==0 || Q_a_max < S) {
                    a_max = a;
                    Q_a_max = S;
                }
            }
            V[s] = Q_a_max;
            dV[s] = pV[s] - V[s];
            pV[s] = V[s];
        }
        Delta = Max(dV) - Min(dV);
    max_iter--;
    
    } while(Delta >= threshold && max_iter > 0);
	
}


/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void OptimisticValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
}
