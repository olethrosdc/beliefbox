// -*- Mode: c++ -*-
// copyright (c) 2005-2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTEXT_BANDIT_GAUSSIAN
#define CONTEXT_BANDIT_GAUSSIAN

#include "MDPModel.h"
#include "NormalDistribution.h"
#include "MeanEstimator.h"
#include "real.h"
#include <vector>

class ContextBanditGaussian : public MDPModel
{
protected:
    std::vector<NormalDistributionUnknownMean> ER;
    int N;
    virtual int getID (int s, int a) const
    {
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }
public:
    ContextBanditGaussian (int n_states, int n_actions, real tau, real mu_0, real tau_0);
    virtual ~ContextBanditGaussian();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual Vector getTransitionProbabilities (int s, int a) const;
    virtual real getRewardDensity(int s, int a, real r) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    virtual void ShowModel() const;
    virtual DiscreteMDP* generate()
    {
        fprintf(stderr, "Not implemented!\n");
        exit(-1);
    }
    virtual const DiscreteMDP* const getMeanMDP() const
    {
        fprintf(stderr, "There is no mean MDP!\n");
        exit(-1);
    }
    //void SetNextReward(int s, int a, real r);
};


#endif

