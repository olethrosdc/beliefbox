// -*- Mode: c++ -*-
// copyright (c) 2005 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: MDPModel.h,v 1.2 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $
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

#define DEBUG_MDP_MODELS


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
  \class MDPModel
  \brief A model of a Markov Decision Process

  Implements algorithms for learning the parameters of an observed
  MDP.  This is an abstract class, but it can only be used for
  discrete spaces.
  
  Since this the _model_ of an MDP, it can actually be used to create
  an MDP.  However, it is _not_ an MDP.

  An implementation of this class is DiscreteMDPCounts

*/
class MDPModel 
{
protected:
    int n_states; ///< number of states (or dimensionality of state space)
    int n_actions; ///< number of actions (or dimensionality of action space)
public:
    MDPModel (int n_states, int n_actions)
    {
        this->n_states = n_states;
        this->n_actions = n_actions;
        mdp_dbg("Creating MDPModel with %d states and %d actions\n",  n_states, n_actions);
    }
    inline virtual int GetNStates()
    {
        return n_states;
    }
    inline virtual int GetNActions()
    {
        return n_actions;
    }
    virtual ~MDPModel()
    {
    }
    virtual void AddTransition(int s, int a, real r, int s2) = 0;
    virtual real GenerateReward (int s, int a) const
    {
        return 0.0;
    }
    virtual int GenerateTransition (int s, int a) const = 0;
    virtual real getTransitionProbability (int s, int a, int s2) const
    {
        return 0.0;
    }
    virtual real getExpectedReward (int s, int a) const
    {
        return 0.0;
    }
    virtual void Reset() = 0;
    virtual DiscreteMDP* CreateMDP();
    virtual DiscreteMDP* generate() = 0;
    virtual const DiscreteMDP* getMeanMDP() const = 0;
    virtual void ShowModel();

};

class GradientDescentMDPModel : public MDPModel
{
protected:
    real alpha; ///< learning rate.
    real** P;
    real* R;
    int N;
    virtual int getID (int s, int a) const
    {
        SMART_ASSERT(s>=0 && s<n_states)(s)(n_states);
        SMART_ASSERT(a>=0 && a<n_actions)(a)(n_actions);
        return s*n_actions + a;
    }
public:
    Distribution* initial_transitions;
    Distribution* initial_rewards; 
    GradientDescentMDPModel (int n_states, int n_actions, Distribution* initial_transitions, Distribution* initial_rewards);
    virtual ~GradientDescentMDPModel();
    virtual void AddTransition(int s, int a, real r, int s2);

    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual real getExpectedReward (int s, int a) const;
    //virtual int getNStates () { return N;}
    virtual void setLearningRate (real learning_rate)
    {
        alpha = learning_rate;
        assert (alpha>=0.0f && alpha <=1.0f);
    }
    virtual void Reset();
};


#endif
