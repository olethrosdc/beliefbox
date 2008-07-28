/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cmath>
#include <limits>
#include "UCB.h"



UCB1Policy::UCB1Policy(int n_actions_) : n_actions(n_actions_)
{
    Reset();
}

UCB1Policy::~UCB1Policy()
{
}

void UCB1Policy::Reset()
{
    times.resize(n_actions);
    Er.resize(n_actions);
    U.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        times[i] = 0;
        Er[i] = 10.0;
        U[i] = INF;
    }
    t = 0;
}

void UCB1Policy::Observe(int a, real r)
{
    real n_j = (real) (++times[a]);
    Er[a] += (r - Er[a]) / n_j;
    t++;
}

int UCB1Policy::SelectAction()
{
    for (int a=0; a<n_actions; a++) {
        if (times[a]==0) {
            return a;
        }
        U[a] = Er[a] + sqrt( 2.0*log((real) t)/ (real) times[a]);
    }
    return ArgMax(n_actions, &U[0]);
}

void UCB1Policy::SetTimes(int a, int times_)
{
    times[a] = times_;
}

void UCB1Policy::SetReward(int a, real reward)
{
    Er[a] = reward;
}

UCBgPolicy::UCBgPolicy(int n_actions_, real gamma_)
    : n_actions(n_actions_), gamma(gamma_)
{
    assert (gamma >= 0 && gamma <= 1);
    Reset();
}

UCBgPolicy::~UCBgPolicy()
{
}

void UCBgPolicy::Reset()
{
    times.resize(n_actions);
    Er.resize(n_actions);
    U.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        times[i] = 0;
        Er[i] = 10.0;
        U[i] = INF;
    }
    t = 0;
}

void UCBgPolicy::Observe(int a, real r)
{
    real n_j = (real) (++times[a]);
    Er[a] += (r - Er[a]) / n_j;
    t++;
}

int UCBgPolicy::SelectAction()
{
    for (int a=0; a<n_actions; a++) {
        if (times[a]==0) {
            return a;
        }
        real k = (real) times[a];
        U[a] = Er[a] + sqrt( 5.0 * log(1.0 / (1 - gamma)) / k);
    }
    return ArgMax(n_actions, &U[0]);
}

void UCBgPolicy::SetTimes(int a, int times_)
{
    times[a] = times_;
}

void UCBgPolicy::SetReward(int a, real reward)
{
    Er[a] = reward;
}


UCB2Policy::UCB2Policy(int n_actions_) : n_actions(n_actions_)
{
    Reset();
}

UCB2Policy::~UCB2Policy()
{
}

void UCB2Policy::Reset()
{
    times.resize(n_actions);
    Er.resize(n_actions);
    U.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        times[i] = 0;
        Er[i] = 0.0;
        U[i] = INF;
    }
    t = 0;
}

void UCB2Policy::Observe(int a, real r)
{
    real n_j = (real) (++times[a]);
    Er[a] += (r - Er[a]) / n_j;
    U[a] = Er[a] + sqrt(2*log((real) t)/n_j);
}

int UCB2Policy::SelectAction()
{
    return ArgMax(n_actions, &U[0]);
}
