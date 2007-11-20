/* -*- Mode: C++; -*- */
/* $Revision */
// copyright (c) 2006,2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTINUOUS_ACTION_POLICY_H
#define CONTINUOUS_ACTION_POLICY_H

#include "AbstractPolicy.h"
#include "ContinuousActionMDPDistribution.h"
#include "Vector.h"
#include <vector>

typedef AbstractPolicy<int, Vector> AbstractContinuousActionPolicy;
typedef Vector Action;

class SamplingContinuousActionPolicy : public AbstractContinuousActionPolicy
{
protected:
    ContinuousActionMDPDistribution& belief;
public:
    SamplingContinuousActionPolicy(ContinuousActionMDPDistribution& belief_);
    virtual ~SamplingContinuousActionPolicy();
    virtual Action SelectAction();
    virtual void Observe (int previous_state,
                          Action action,
                          real r,
                          int next_state);
    virtual void Reset();

};

class SphereContinuousActionPolicy : public AbstractContinuousActionPolicy
{
protected:
    int d;
public:
    SphereContinuousActionPolicy(int d);
    virtual ~SphereContinuousActionPolicy();
    virtual Action SelectAction();
    virtual void Observe (int previous_state,
                          Action action,
                          real r,
                          int next_state);
    virtual void Reset();
};

#endif
