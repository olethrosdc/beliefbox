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

#ifndef DISCRETE_MDP_AGGREGATE
#define DISCRETE_MDP_AGGREGATE

#include "DiscreteMDPCounts.h"
#include <set>
#include <vector>

/// A container for states
class DiscreteStateAggregate
{
public:
    std::set<int> S;
    void add(int state)
    {
        S.insert(state);
    }
    bool contains(int state)
    {
        return (S.find(state) != S.end());
    }
};


/// USes DiscreteStateAggregate to form an MDP
class DiscreteMDPAggregate : public DiscreteMDPCounts
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
public:
    DiscreteMDPAggregate (int n_aggregated_states, int n_states, int n_actions, int init_transition_count=0, int init_reward_count = 0, real init_reward = 0.0);
    virtual ~DiscreteMDPAggregate();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual Vector getTransitionProbabilities (int s, int a) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    void SetNextReward(int s, int a, real r);
};


#endif
