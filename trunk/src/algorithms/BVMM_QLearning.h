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
#ifndef BVMM_QLEARNING_H
#define BVMM_QLEARNING_H

#include "real.h"
#include "debug.h"
#include "OnlineAlgorithm.h"
#include <cstdio>

/** Q-Learning on a Bayesian Variable Order Markov Model

    From the paper:
    "Context model inference for large or partially observable MDPs",
    C. Dimitrakakis, ICML workshop on RL and search, 2010.
 */
template <class T>
class BVMM_QLearning : public OnlineAlgorithm<int, int>
{
protected:
    int n_actions; ///< the number of actions
    int n_obs; ///< the number of distinct observations
    T tree; ///< the context tree
    int n_steps; ///< number of steps passed
    int n_episode_steps; ///< number of steps passed in an episode
    int action; ///< the current action
    int current_obs; ///< the current observation
    real current_reward; ///< the current reward
    real epsilon; ///< the randomness of the action selection
    real alpha; ///< the step-size
    real gamma; ///< the eligibility trace
public:
    BVMM_QLearning(int n_actions_, int n_obs_, int depth, real epsilon_, real alpha_, real gamma_)
        : n_actions(n_actions_),
          n_obs(n_obs_),
          tree(n_obs * n_actions, n_obs, n_actions, n_obs, depth),
          n_steps(0),
          n_episode_steps(0),
          current_obs(0),
          epsilon(epsilon_),
          alpha(alpha_),
          gamma(gamma_)
    {        
    }
    // --------- OnlineAlgorithm functions -----------
    virtual ~BVMM_QLearning()
    {
    }
    
    /// call this at the end of an episode.
    virtual void Reset()
    {
        n_episode_steps = 0;
        action = 0;
    }
    /// Partial observation
    virtual real Observe (real reward, int next_state, int next_action)
    {
        _Observe(next_action, next_state, reward);
        return _QLearning(alpha, gamma);
    }
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, int next_state) 
    {
        assert(next_state >= 0 && next_state < n_obs);

        if (n_episode_steps) {
            _Observe(action, next_state, reward);
        } else {
            _Observe(next_state);
        }
        n_episode_steps++;
        n_steps++;
        //_QLearning(alpha, gamma);
        real threshold = epsilon / (1.0 + sqrt(n_steps));
        if (urandom() < threshold) {
            action =  rand()%n_actions;
            return action;
        }
        real Q_max = QValue(0);
        int arg_max = 0;
        for (int i=0; i<n_actions; ++i) {
            real Q = QValue(i);
            //printf ("%f ", Q);
            if (Q > Q_max) {
                Q_max = Q;
                arg_max = i;
            }
        }
        //printf(" # Q\n");
        action = arg_max;
        return action;
    }

    ///
    virtual real getValue (int state, int action) 
    {
        return QValue(action);
    }

    
    // --------- Factored Algorithm hooks -------------
    /* Training and generation */
    /// Observe the (first?) observation.
    virtual real _Observe (int prd) 
    {
        current_obs = prd;
        return 1.0 / (real) n_obs;
    }
    
    /// Observe current action and next observation
    real _Observe (int act, int prd, real reward) 
    {
        int x = act * n_obs + current_obs;
        current_obs = prd;
        current_reward = reward;
        //printf ("tree: %d %d %f\n", x, prd, reward);
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
    virtual real _QLearning(real step_size, real gamma)
    {
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
};

#endif
