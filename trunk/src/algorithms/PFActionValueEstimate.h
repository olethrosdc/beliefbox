/* -*- Mode: C++; -*- */
/* VER: $Id: PFActionValueEstimate.h,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PF_ACTION_VALUE_ESTIMATE
#define PF_ACTION_VALUE_ESTIMATE

#include "ActionValueEstimate.h"
#include "Distribution.h"

/** A population estimate of actions
 */
class PFActionValueEstimate : public ActionValueEstimate
{
public:
    std::vector<BernoulliGridParticleFilter> q; ///< The estimates
    std::vector<real> s; ///< The samples
    int n_actions;
    UniformDistribution* transitions;
    BernoulliDistribution* observations;
    UniformDistribution* prior;
    PFActionValueEstimate(int n_actions, int n_members);
    virtual void Reset();
    virtual ~PFActionValueEstimate();
    virtual void Observe (int a, real r);
    virtual real GetMean (int a);
    virtual real GetVar (int a);
    virtual int Sample();
    virtual int GetMax();
    virtual int GetSecondMax();
    virtual real Sample(int a);
    virtual real GetMemberValue(int a, int i);
    virtual real GetProbability(int i, int j, real delta);
};

#endif





