/* -*- Mode: C++; -*- */
/* VER: $Id: SimpleStateActionEstimator.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SimpleStateActionEstimator.h"
#include "real.h"
#include "SmartAssert.h"

#include <cstdlib>
#include <cstdio>


SimpleStateActionEstimator::SimpleStateActionEstimator(int n_states, int n_actions, real gamma, real init_val) : StateActionEstimator(n_states, n_actions)
{
    this->init_val = init_val;

    Q_data = new real [n_actions*n_states];
    Q = new real* [n_states];
    for (int i=0; i<n_states; i++) {
        Q[i] = &Q_data[i*n_actions];
    }
    Reset();
}

SimpleStateActionEstimator::~SimpleStateActionEstimator()
{
    delete [] Q;
    delete [] Q_data;
}

void SimpleStateActionEstimator::Reset()
{
    for (int i=0; i<n_states; i++) {
        for (int j=0; j<n_states; j++) {
            Q[i][j]=init_val;
        }
    }
}

real SimpleStateActionEstimator::TransProbability(int s_p, int a, int s)
{
    fprintf (stderr, "SimpleStateActionEstimator cannot estimate transition probability\n");
    exit(-1);
}

int SimpleStateActionEstimator::SampleTransition(int s, int a)
{
    fprintf (stderr, "SimpleStateActionEstimator cannot estimate transition probability\n");
    exit(-1);
}

real SimpleStateActionEstimator::GetTransitionPrior(int s_p, int a, int s, real x)
{
    fprintf (stderr, "SimpleStateActionEstimator cannot estimate transition probability\n");
    exit(-1);
}

real SimpleStateActionEstimator::SampleTransitionPrior(int s_p, int a, int s)
{
    fprintf (stderr, "SimpleStateActionEstimator cannot estimate transition probability\n");
    exit(-1);
}


//-------- MaxStateActionEstimator -------------//
/** Estimator of maximum return following a state-action pair
 *
 *  This is very useful for Q-learning.
 */
MaxStateActionEstimator::MaxStateActionEstimator(int n_states, int n_actions, real alpha, real gamma, real init_val) : SimpleStateActionEstimator(n_states, n_actions, gamma, init_val)
{
    this->alpha = alpha;
}

/// Destructor
MaxStateActionEstimator::~MaxStateActionEstimator()
{
}


/** Obesrve a (s,a,r,s') tuplet
 *
 * Update max estimate of s_p,a_p, ignoring next state.
 */
void MaxStateActionEstimator::Observe(int s_p, int a_p, real r, int s, int a)
{
    //DISABLED_ASSERT(s>=0);
    //DISABLED_ASSERT(a>=0);
    real Er = MaxActionValue(s);
    real TD = r + gamma*Er - Q[s_p][a_p];
    Q[s_p][a_p] += alpha * TD;
}


//-------- SarsaStateActionEstimator -------------//
/** Estimator of expected return.
 *
 *  This is very useful for SARSA.
 */
SarsaStateActionEstimator::SarsaStateActionEstimator(int n_states, int n_actions, real alpha, real gamma, real init_val) : SimpleStateActionEstimator(n_states, n_actions, gamma, init_val)
{
    this->alpha = alpha;
}

/** Obesrve a (s,a,r,s') tuplet
 *
 * Update estimate of s_p,a_p under current policy
 */
void SarsaStateActionEstimator::Observe(int s_p, int a_p, real r, int s, int a)
{
    //DISABLED_ASSERT(s>=0);
    //DISABLED_ASSERT(a>=0);
    real Er = Q[s][a];
    real TD = r + gamma*Er - Q[s_p][a_p];
    Q[s_p][a_p] += alpha * TD;
}
