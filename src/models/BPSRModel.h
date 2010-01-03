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

#ifndef BPSR_MODEL_H
#define BPSR_MODEL_H

#include <cmath>
#include <cassert>
#include <vector>
#include "real.h"
#include "SmartAssert.h"
#include "Distribution.h"
#include "DiscreteMDP.h"
#include "BayesianPredictiveStateRepresentation.h"
#include "DiscreteVariable.h"

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
    BayesianPredictiveStateRepresentation* bpsr;
    DiscreteVector* Z;
    std::vector<int> getIndexVector (int x, real r) 
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
        std::vector<int> v(2);
        v[0] = x;
        v[1] = r_i;
        return v;
    }
public:
    BPSRModel  (int n_obs_, int n_actions_, std::vector<real> rewards_, int tree_depth);

    virtual ~BPSRModel();
    virtual void Observe(int a, int x, real r);
    virtual real getTransitionProbability(std::vector<int> history, int a, int x) const;
    virtual real getTransitionProbability(int a, int x) const;
    virtual real getExpectedReward (std::vector<int> history) const;
    virtual void Reset();
};


#endif
