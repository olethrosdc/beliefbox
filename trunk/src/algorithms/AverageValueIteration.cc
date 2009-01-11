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

#include "AverageValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include "Random.h"
#include <cmath>
#include <cassert>

AverageValueIteration::AverageValueIteration(const DiscreteMDP* mdp, real baseline)
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

void AverageValueIteration::Reset()
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
    p_b.resize(n_states);
    real sum = 0.0;
    for (int s=0; s<n_states; s++) {
        p_b[s] = urandom();
        sum += p_b[s];
    }
    for (int s=0; s<n_states; s++) {
        p_b[s] /= sum;
    }
}

AverageValueIteration::~AverageValueIteration()
{
}

void AverageValueIteration::ComputeStateValues(real threshold, int max_iter)
{
    std::vector<real> p_tmp(n_states);
    std::vector<int> a_max(n_states);
    int n_policy_changes = 0;
    do {
        Delta = 0.0;
#if 0
        baseline = 0.0;
        for (int s=0; s<n_states; s++) {
            baseline += p_b[s] * V[s];
        }
#else
        baseline = 0.0;
        for (int s=0; s<n_states; s++) {
            if (baseline < V[s]) {
                baseline = V[s];
            }
        }
#endif
        //baseline = = real(n_states);
        n_policy_changes = 0;
        for (int s=0; s<n_states; s++) {
            real Q_a_max = -RAND_MAX;
            int c_a_max = 0;
            for (int a=0; a<n_actions; a++) {
                real S = 0.0;
                DiscreteStateSet next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                    real R = mdp->getExpectedReward(s, a) + V[s2] - baseline;
                    S += P * R;
                }
                if (a==0 || Q_a_max < S) {
                    c_a_max = a;
                    Q_a_max = S;
                }
            }
            if (c_a_max != a_max[s]) {
                a_max[s] = c_a_max;
                n_policy_changes++;
            }
            V[s] = Q_a_max;
            dV[s] = pV[s] - V[s];
            pV[s] = V[s];
        }

        for (int s=0; s<n_states; s++) {
            p_tmp[s] = 0.0;
        }
        for (int s=0; s<n_states; s++) {
            // calculate new p_b
            DiscreteStateSet next = mdp->getNextStates(s, a_max[s]);
            for (DiscreteStateSet::iterator i=next.begin();
                 i!=next.end();
                 ++i) {
                int s2 = *i;
                real P = mdp->getTransitionProbability(s, a_max[s], s2);
                p_tmp[s2] += 0.01 +  p_b[s] * P;
            }
        }
#if 1
        real sum = 0.0;
        for (int s=0; s<n_states; s++) {
            sum += p_tmp[s];
        }
        //printf ("sum: %f\n", sum);
        for (int s=0; s<n_states; s++) {
            p_b[s] = p_tmp[s] / sum;
            //p_b[s] = 1/((real) n_states);//xp_tmp[s] / sum;
        }
#else
        for (int s=0; s<n_states; s++) {
            p_b[s] = p_tmp[s];
        }
#endif
        Delta = Max(dV) - Min(dV);
        max_iter--;

    } while(Delta >= threshold && max_iter > 0 && n_policy_changes > 0);
	if (n_policy_changes == 0) {
        printf ("Policy did not change : iters left = %d\n", max_iter);
    }
}


/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void AverageValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
    int N = n_states * n_actions;

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            dQ[s][a] = 0.0;
        }
    }
    std::vector<real> p_b(n_states);
    real sum = 0.0;
    for (int s=0; s<n_states; s++) {
        p_b[s] = urandom();
        sum += p_b[s];
    }
    for (int s=0; s<n_states; s++) {
        p_b[s] /= sum;
    }
    do {
        Delta = 0.0;
        baseline = 0.0;
        for (int s=0; s<n_states; s++) {
            baseline += p_b[s] * V[s];
        }
        baseline = 0.0;///= real(n_states);
        for (int s0=0; s0<n_states; s0++) {
            int s = s0;
            for (int a=0; a<n_actions; a++) {
                real sum = 0.0;

                DiscreteStateSet next = mdp->getNextStates(s, a);
                for (DiscreteStateSet::iterator i=next.begin();
                     i!=next.end();
                     ++i) {
                    int s2 = *i;
                    real P = mdp->getTransitionProbability(s, a, s2);
                        //if (P > 0) {
                        real R = mdp->getExpectedReward(s, a) - baseline;
                        real Q_a_max = Max(n_actions, Q[s2]);
                        sum += P*(R + gamma*Q_a_max);
                            //}

                }
                Q[s][a] = sum;
                dQ[s][a] = pQ[s][a] - sum;
                pQ[s][a] = sum;
                //baseline += 0.1 *(V[s] - baseline);
            }
        }
        
        Delta = Max(N, &dQ_data[0]) - Min(N, &dQ_data[0]);			
        max_iter--;
    } while(Delta >= threshold && max_iter > 0);
}

