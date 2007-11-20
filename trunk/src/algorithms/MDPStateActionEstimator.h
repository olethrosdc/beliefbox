/* -*- Mode: C++; -*- */
/* VER: $Id: MDPStateActionEstimator.h,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MDP_STATE_ACTION_ESTIMATOR_H
#define MDP_STATE_ACTION_ESTIMATOR_H

#include "MDPModel.h"
#include "StateActionEstimator.h"
#include "MathFunctions.h"
#include "real.h"


/** A point estimate for actions
 *   
 *  This is a base class for all estimators, such as E[R|pi] or
 *  for off-policy estimates such as E[R|pi^*].
 */
class MDPStateActionEstimator : public StateActionEstimator
{
 public:
    MDPModel* model;
    MDPStateActionEstimator(MDPModel* model, real gamma, real init_val = 0.0f);
    virtual ~MDPStateActionEstimator();
    virtual void Reset();
    virtual void Observe(int s_p, int a_p, real r, int s=-1, int a=-1) = 0;
    virtual real TransProbability(int s_p, int a_p, int s);
    virtual int SampleTransition(int s, int a);
    virtual real GetExpectedReward(int s, int a);
    virtual real SampleReward(int s, int a);

};

#if 0
/** A prior estimate for actions
 *   
 *  This is a base class for all estimators, such as E[R|pi] or
 *  for off-policy estimates such as E[R|pi^*].
 */
class PriorMDPEstimator : public MDPStateActionEstimator
{
 public:
    PriorMDPEstimator();
    virtual ~PriorMDPEstimator();
};
#endif

#endif
