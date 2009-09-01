/* -*- Mode: C++; -*- */
/* VER: $Id: SimpleStateActionEstimator.h,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef SIMPLE_STATE_ACTION_ESTIMATOR_H
#define SIMPLE_STATE_ACTION_ESTIMATOR_H

#include "StateActionEstimator.h"
#include "MathFunctions.h"
#include "real.h"

/** A point estimate for actions
 *   
 *  This is a base class for all estimators, such as E[R|pi] or
 *  for off-policy estimates such as E[R|pi^*].
 */
class SimpleStateActionEstimator : public StateActionEstimator
{
 public:
    real* Q_data; ///< estimated return
    real** Q;  ///< estimated return
    real init_val; ///< initial value for estimates
    real alpha;  ///< learning rate
    real gamma; ///< gamma
    SimpleStateActionEstimator(int n_states, int actions, real gamma, real init_val = 0.0f);
    virtual ~SimpleStateActionEstimator();
    virtual void Reset();
    virtual void Observe(int s_p, int a_p, real r, int s=-1, int a=-1) = 0;
    virtual real TransProbability(int s_p, int a, int s);
    virtual int SampleTransition(int s, int a) ;
    virtual real GetTransitionPrior(int s_p, int a, int s, real x);
    virtual real SampleTransitionPrior(int s_p, int a, int s);

};

class MaxStateActionEstimator : public SimpleStateActionEstimator
{
 public:
    MaxStateActionEstimator(int n_states, int actions, real alpha, real gamma, real init_val = 0.0f);
    virtual ~MaxStateActionEstimator();
    virtual void Observe(int s_p, int a_p, real r, int s=-1, int a=-1);
    real MaxActionValue(int state)
    {
        return Max(n_actions, Q[state]);
    }
};

class SarsaStateActionEstimator : public SimpleStateActionEstimator
{
 public:
    SarsaStateActionEstimator(int n_actions, int n_states, real alpha, real gamma, real init_val = 0.0f);
    virtual ~SarsaStateActionEstimator();
    virtual void Observe(int s_p, int a_p, real r, int s=-1, int a=-1);
};



#endif
