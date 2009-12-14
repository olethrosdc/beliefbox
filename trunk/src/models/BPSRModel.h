// -*- Mode: c++ -*-
// copyright (c) 2005 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: BPSRModel.h,v 1.2 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $
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
   \ingroup MachineLearning
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
    BayesianPredictiveStateRepresentation bpsr;
public:
    BPSRModel  (int n_obs, int n_actions)
    {
        this->n_states = n_states;
        this->n_actions = n_actions;
        mdp_dbg("Creating BPSRModel with %d states and %d actions\n",  n_states, n_actions);
    }
    virtual ~BPSRModel()
    {
    }
    virtual void AddTransition(int a, int r, int s) = 0;
    virtual real GenerateReward (int s, int a) const
    {
        return 0.0;
    }
    virtual int GenerateTransition (int s, int a) const = 0;
    virtual real getTransitionProbability (int s, int a, int s2) const
    {
        return 0.0;
    }
    virtual real getTransitionProbability(std::vector<int> history, int a, int x);
    virtual real getExpectedReward (int s, int a) const
    {
        return 0.0;
    }
    virtual void Reset() = 0;
    virtual DiscreteMDP* CreateMDP();
    virtual const DiscreteMDP* getMeanMDP() const = 0;
    virtual void ShowModel();

};


#endif
