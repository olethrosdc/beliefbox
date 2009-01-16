// -*- Mode: c++ -*-
// copyright (c) 2005-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_MDP_H
#define DISCRETE_MDP_H

#include "SmartAssert.h"
#include "MDP.h"
#include "Distribution.h"
#include <vector>
#include <set>

typedef std::set<int> DiscreteStateSet;
typedef std::set<int>& DiscreteStateSetRef;

template<>
class MDP<int, int> {
protected:
    int state;
    int n_states; ///< number of states (or dimensionality of state space)
    int n_actions; ///< number of actions (or dimensionality of action space)
    std::vector<real*> P; ///< transition distribution
    std::vector<real> P_data; ///< transition distribution data
    std::vector<Distribution*> R; ///< reward distribution
    std::vector<real> ER; ///< expected reward
    std::vector<DiscreteStateSet> next_states;
    int N;
    inline int getID (int s, int a) const
    {
#ifndef NDEBUG
        SMART_ASSERT(s>=0 && s<n_states)(s)(n_states);
        SMART_ASSERT(a>=0 && a<n_actions)(a)(n_actions);
#endif
        return s*n_actions + a;
    }

public:

    MDP<int, int>(int n_states, int n_actions, real** initial_transitions, Distribution** initial_rewards);
    MDP<int,int> (const MDP<int,int>& mdp);

    inline int GetNStates() const
    {
        return n_states;
    }
    inline int GetNActions() const
    {
        return n_actions;
    }
    virtual ~MDP<int, int>();
    virtual void ShowModel() const;
    virtual void dotModel(FILE* fout) const;
    real generateReward (int s, int a) const;
    int generateState (int s, int a) const;
    inline real getTransitionProbability (int s, int a, int s2) const
    {
        int ID = getID (s, a);                
        assert (s2>=0 && s2<n_states);
        return P[ID][s2];
    }
    inline real getExpectedReward (int s, int a) const
    {
        int ID = getID (s, a);
        return ER[ID];
    }
    inline void setTransitionProbability(int s, int a, int s2, real p)
    {
        int ID = getID (s, a);
        real* Ps=P[ID];
        assert(s2>=0 && s2<n_states);
        Ps[s2] = p;
        DiscreteStateSet& next = next_states[ID];
        if (p==0) {
            next.erase(s2);
        } else {
            next.insert(s2);
        }
    }
    inline void setRewardDistribution(int s, int a, Distribution* reward)
    {   
        int ID = getID (s, a);
        R[ID] = reward;
        ER[ID] = reward->getMean();
    }
    inline DiscreteStateSet getNextStates(int s, int a) const
    {
        int ID = getID (s,a);
        return next_states[ID];
    }
    void AperiodicityTransform(real tau);
    void Check() const;
    real CalculateDiameter() const;
};

typedef MDP<int, int> DiscreteMDP;

#endif
