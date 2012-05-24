// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTEXT_BANDIT_COLLECTION
#define CONTEXT_BANDIT_COLLECTION

#include "ContextBanditGaussian.h"
#include "ContextBanditAggregate.h"
#include "Gridworld.h"
#include <set>
#include <vector>

/// blah
class ContextBanditCollection : public MDPModel
{
protected:
    std::vector<ContextBanditGaussian*> A;
    std::vector<real> P;
public:
    ContextBanditCollection(int n_aggregates, int n_states, int n_actions, real tau, real mu_0, real tau_0);
    virtual ~ContextBanditCollection();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual real getRewardDensity (int s, int a, real r) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    virtual DiscreteMDP* generate() const
    {
        fprintf(stderr, "Not implemented!\n");
        exit(-1);
    }
    virtual const DiscreteMDP* const getMeanMDP() const
    {
        fprintf(stderr, "There is no mean MDP!\n");
        exit(-1);
    }

};



#endif
