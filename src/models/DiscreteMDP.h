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

template<>
class MDP<int, int> {
protected:
	int state;
    int n_states; ///< number of states (or dimensionality of state space)
    int n_actions; ///< number of actions (or dimensionality of action space)
    real** P; ///< transition distribution
    real* P_data; ///< transition distribution data
    Distribution** R; ///< reward distribution
    real* ER; ///< expected reward
    int N;
    int getID (int s, int a)
    {
        SMART_ASSERT(s>=0 && s<n_states)(s)(n_states);
        SMART_ASSERT(a>=0 && a<n_actions)(a)(n_actions);
        return s*n_actions + a;
    }
public:
    MDP<int, int> (int n_states, int n_actions, real** initial_transitions, Distribution** initial_rewards);

    int GetNStates()
    {
        return n_states;
    }
    int GetNActions()
    {
        return n_actions;
    }
    virtual ~MDP<int, int>();
    virtual void ShowModel();
    virtual real generateReward (int s, int a);
    virtual int generateState (int s, int a);
    virtual real getTransitionProbability (int s, int a, int s2);
  //	virtual real getRewardProbability (int s, int a);
    virtual real getExpectedReward (int s, int a);
    virtual void setTransitionProbability(int s, int a, int s2, real p);
    virtual void setRewardDistribution(int s, int a, Distribution* reward);
    virtual void Check();
};

typedef MDP<int, int> DiscreteMDP;

#endif
