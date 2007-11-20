/* -*- Mode: C++; -*- */
/* VER: $Id: PFActionValueEstimate.c,v 1.2 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "PFActionValueEstimate.h"

PFActionValueEstimate::PFActionValueEstimate(int n_actions, int n_members)
{
    transitions = new UniformDistribution (-0.1f, 0.1f);
    prior = new UniformDistribution(0.0, 1.0f);
    observations = new BernoulliDistribution;
    this->n_actions = n_actions;
    q.resize(n_actions);
    s.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        q[i].Init(n_members, prior, transitions, observations);
        q[i].Reset();
    }
}



void PFActionValueEstimate::Reset()
{
    for (int i=0; i<n_actions; i++) {
        q[i].Reset();
    }
}

/// Destroy
PFActionValueEstimate::~PFActionValueEstimate()
{
    delete prior;
    delete observations;
    delete transitions;
}

/// Update estimates when observing reward r after action a is taken
void PFActionValueEstimate::Observe (int a, real r)
{
    q[a].Observe(r);
}

/// Get the mean reward of action a
real PFActionValueEstimate::GetMean (int a)
{
    return q[a].GetMean();
}

/// Get the mean reward of action a
real PFActionValueEstimate::GetVar (int a)
{
    return q[a].GetVar();
}

/// Sample an action i with probability p(q[i] > q[j])
int PFActionValueEstimate::Sample()
{
    for (int i=0; i<n_actions; i++) {
        s[i] = q[i].Sample();
    }
    int arg_max = 0;
    real max = s[0];
    for (int i=1; i<n_actions; i++) {
        if (max < s[i]) {
            max = s[i];
            arg_max = i;
        }
    }
    return arg_max;
}
/// Sample the action with maximum mean
int PFActionValueEstimate::GetMax()
{
    int arg_max = 0;
    real max = q[0].GetMean();
    for (int i=1; i<n_actions; i++) {
        real m = q[i].GetMean();
        if (max < m) {
            max = m;
            arg_max = i;
        }
    }
    return arg_max;
}
int PFActionValueEstimate::GetSecondMax()
{
    int arg_max = GetMax();
    int arg_max2 = 0;
    real max;
    if (arg_max==0) {
        max = q[1].GetMean();
    } else {
        max = q[0].GetMean();
    }
    for (int i=0; i<n_actions; i++) {
        if (i!=arg_max){
            real m = q[i].GetMean();
            if (max < m) {
                max = m;
                arg_max2 = i;
            }
        }
    }
    return arg_max2;
}

real PFActionValueEstimate::Sample(int a)
{
    return q[a].Sample();
}


real PFActionValueEstimate::GetMemberValue(int a, int i)
{
    return q[a].y[i];
}

/// Get a simple estimate of P(q_i - q_j > delta)
real PFActionValueEstimate::GetProbability(int i, int j, real delta)
{
    real P = 0.0f;
    int N = q[i].N;
    int M = q[j].N;
#if 1
    real S = 0.0f;
    for (int n=0; n<N; n++) {
        for (int m=0; m<M; m++) {
            real pair_density = q[i].w[n] * q[j].w[m];
            S += pair_density;
            if (q[i].y[n] - q[j].y[m] > delta) {
                P += pair_density;
            }
        }
    }
    return P / S;
#else
    // pretend to be sampling
    real S=0.0f;
    for (int n=0; n<N; n++) {
        int m = rand()%M;
        S += q[i].w[n];
        if (q[i].y[n] - q[j].y[n] > delta) {
            P += q[i].w[n] *  q[j].w[m];
        }
    }
    return P/S;
#endif
}








