/* -*- Mode: c++;  -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef FACTORED_PREDICTOR_RL_H
#define FACTORED_PREDICTOR_RL_H

#include "real.h"
#include "debug.h"
#include <cstdio>


/// Abstract class for prediction with actios
class FactoredPredictorRL
{
public:
    virtual ~FactoredPredictorRL()
    {}
    
    /* Training and generation */
    virtual real Observe (int prd) = 0;
    virtual real Observe (int act, int prd, real r) = 0;
    virtual real ObservationProbability (int act, int x) = 0;
    virtual real QValue (int act) = 0;
    virtual real QLearning (real step_size, real gamma ) = 0;
    virtual real Sarsa (real step_size, real gamma ) = 0;
    //virtual real ObservationProbability (int x) = 0;
    virtual void Reset() = 0;
    
}; 

template <class T>
class TFactoredPredictorRL : public FactoredPredictorRL
{
protected:
    int n_actions; ///< the number of actions
    int n_obs; ///< the number of distinct observations
    T tree; ///< the context tree
    int current_obs; ///< the current observation
    real current_reward; ///< the current reward
public:
    TFactoredPredictorRL(int n_actions_, int n_obs_, int depth)
        : n_actions(n_actions_),
          n_obs(n_obs_),
          tree(n_obs * n_actions, n_obs, n_actions, n_obs, depth),
          current_obs(0)
    {        
    }

    virtual ~TFactoredPredictorRL()
    {
    }
    /* Training and generation */
    /// Observe the (first?) observation.
    virtual real Observe (int prd) 
    {
        current_obs = prd;
        return 1.0 / (real) n_obs;
    }
    /// Observe current action and next observation
    virtual real Observe (int act, int prd, real reward) 
    {
        int x = act * n_obs + current_obs;
        current_obs = prd;
        current_reward = reward;
        //printf ("%d %d %f\n", x, prd, reward);
        return tree.Observe(x, prd, reward);
    }
    
    /// Observe current action and next observation
    virtual real QValue (int act) 
    {
        //Serror("Not implemented\n");
        int x = act * n_obs + current_obs;
        return tree.QValue(x);
    }


    /// Do q-learning, starting with next observation
    virtual real QLearning(real step_size, real gamma)
    {
        //Serror("Not implemented\n");
        return tree.QLearning(step_size, gamma, current_obs, current_reward);
    }

    /// Do q-learning, starting with next observation
    virtual real Sarsa(real step_size, real gamma)
    {
        //Serror("Not implemented\n");
        return tree.Sarsa(step_size, gamma, current_obs, current_reward);
    }

    virtual real ObservationProbability (int act, int x) 
    {
        Serror("Not implemented\n");
        return -1;
    }

    virtual void Reset()
    {
    }
};

#endif
