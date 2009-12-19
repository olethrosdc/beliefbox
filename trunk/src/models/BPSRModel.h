// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MDP_MODEL_H_
#define MDP_MODEL_H_

#include <cmath>
#include <cassert>
#include <vector>
#include "real.h"
#include "SmartAssert.h"
#include "Distribution.h"
#include "DiscreteMDP.h"
#include "BayesianPredictiveStateRepresentation.h"

#undef DEBUG_MDP_MODELS


#ifdef DEBUG_MDP_MODELS
#define mdp_dbg  printf("%s:%d: %s(): ", __FILE__, __LINE__, __FUNCTION__); printf
#else
//#define mdp_dbg (void)
inline void mdp_dbg(...)
{
    
}
#endif


/**
   \ingroup ReinforcementLearning
*/

/*! 
  \class BPSRModel
  \brief A model of a Partially observable Markov decision process



*/
class BPSRModel
{
protected:
    int n_obs;
    int n_actions;
    std::vector<real> rewards;
    int n_rewards;
    BayesianPredictiveStateRepresentation bpsr;
    DiscreteVector* Z;
    int getIndex (int a, int x, int r) 
    {
            // find closest reward index
        int r_i = 0;
        real r_d = fabs(rewards[0] - r);
        for (uint i=1; i<rewards.size(); ++i) {
            real d = fabs(rewards[i] - r);
            if (d < r_d) {
                r_i = i;
                r_d = d;
            }
        }

            // set up vector
        std::vector<int> v(3);
        v[0] = a;
        v[1] = x;
        v[2] = r_d;
    }
public:
    BPSRModel  (int n_obs_, int n_actions_, std::vector<real> rewards_)
        : n_obs(n_obs_),
          n_actions(n_actions_),
          rewards(rewards_)
    {
        mdp_dbg("Creating BPSRModel with %d states, %d actions and %d rewards\n",  n_states, n_actions, rewards.size());
        n_rewards = size(rewards);
        std::vector<int> sizes(3);
        sizes(0) = n_obs;
        sizes(1) = n_actions;
        sizes(2) = n_rewards;
        Z = new DiscreteVector(x);
        Z->size();
    }

    virtual ~BPSRModel()
    {
    }

    virtual void AddTransition(int a, int r, int x) 
    {
        
    }
        /// P(x_{t+1} | a_{t+1} , z^t)
    virtual real getTransitionProbability(std::vector<int> history, int a, int x) const;
    virtual real getExpectedReward (std::vector<int> history) const;
    virtual void Reset();
};


#endif
