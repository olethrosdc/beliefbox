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

#ifndef CONTEXT_BANDIT_AGGREGATE
#define CONTEXT_BANDIT_AGGREGATE

#include "ContextBanditGaussian.h"
#include "DiscreteStateAggregate.h"

#include <vector>




/// USes DiscreteStateAggregate to form an MDP
class ContextBanditAggregate : public ContextBanditGaussian
{
protected:
    std::vector<int> state_map; ///< each entry tells you which set the state belongs in
    std::vector<DiscreteStateAggregate> X; ///< each X contains the set of states that are aggregated
    int n_aggregated_states; // the original number of states
    /// Return the Aggregate state in which s belongs
    virtual int getID (int x, int a) const
    {
        int s = Aggregate(x);
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }
    int Aggregate(int s) const
    {
        return state_map[s];
    }
    void BuildRandomAggregate();
    void BuildFactoredAggregate(int n_factors, int n_values);
public:
    ContextBanditAggregate (int n_aggregated_states, int n_states, int n_actions, real tau, real mu_0, real tau_0);
    ContextBanditAggregate (bool blah, int n_factors, int n_values, int n_aggregated_states, int n_states, int n_actions, real tau, real mu_0, real tau_0);

    virtual ~ContextBanditAggregate();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getRewardDensity (int s, int a, real r) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual Vector getTransitionProbabilities (int s, int a) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    virtual DiscreteMDP* generate() 
    {
        fprintf(stderr, "Not implemented!\n");
        exit(-1);
    }
    virtual const DiscreteMDP* getMeanMDP() const
    {
        fprintf(stderr, "There is no mean MDP!\n");
        exit(-1);
    }
    //void SetNextReward(int s, int a, real r);
};
 

#endif
