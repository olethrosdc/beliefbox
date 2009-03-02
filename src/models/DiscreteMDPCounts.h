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

#ifndef DISCRETE_MDP_COUNTS
#define DISCRETE_MDP_COUNTS

#include "MDPModel.h"
#include "Dirichlet.h"
#include "MeanEstimator.h"
#include "real.h"
#include <vector>

class DiscreteMDPCounts : public MDPModel
{
protected:
    std::vector<DirichletDistribution> P;
    std::vector<MeanEstimator> ER;
    int N;
    virtual int getID (int s, int a) const
    {
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }
public:
    DiscreteMDPCounts (int n_states, int n_actions, int init_transition_count=1, int init_reward_count = 0, real init_reward = 0.0);
    virtual ~DiscreteMDPCounts();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual Vector getTransitionProbabilities (int s, int a) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    virtual void ShowModel() const;
    //void SetNextReward(int s, int a, real r);
};


#endif

