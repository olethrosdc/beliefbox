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
#include "Vector.h"
#include <cstdio>
#include <cassert>

/// Abstract class for prediction with actios
template <class X, class A>
class FactoredPredictorRL
{
public:
    virtual ~FactoredPredictorRL()
    {}
    
    /* Training and generation */
	/// Observe a new observation X
    virtual real Observe (const X& obs) = 0;
	/// Observe the current action \c act and the next observation \c obs and reward \c r.
    virtual real Observe (const A& act, const X& obs, real r) = 0;
	/// Return the probability (density) of observing \c x aftering taking action \c act.
    virtual real ObservationProbability (const A& act, const X& x) = 0;
	/// Return the expected utility of taking action \c act.
    virtual real QValue (const A& act) = 0;
	/// Perform a QLearning step
    virtual real QLearning (real step_size, real gamma ) = 0;
	/// Perform a Sarsa step
    virtual real Sarsa (real step_size, real gamma, real epsilon = 0.01 ) = 0;
    //virtual real ObservationProbability (int x) = 0;
    virtual void Reset() = 0;
    virtual void Show()
    {
        printf("FactoredPredictorRL: Nothing to show\n");
    }
    
}; 

/// Abstract class for prediction with actions
template <class T>
class TFactoredPredictorRL : public FactoredPredictorRL<int, int>
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
    virtual real Observe (const int& obs) 
    {
        current_obs = obs;
        return 1.0 / (real) n_obs;
    }
    /// Observe current action and next observation and reward.
    ///
    /// As a side-effect, the current observation changes.
    virtual real Observe (const int& act, const int& obs, real reward) 
    {
        int x = act * n_obs + current_obs;
        current_obs = obs;
        current_reward = reward;
        //printf ("%d %d %f\n", x, obs, reward);
        return tree.Observe(x, obs, reward);
    }
    
    /// Return expected utility of action \c act.
    virtual real QValue (int act) 
    {
        assert (act >= 0 && act < n_actions);
        assert (current_obs >= 0 && current_obs < n_obs);
        //Serror("Not implemented\n");
        int x = act * n_obs + current_obs;
        return tree.QValue(x);
    }

    /// Do Q-learning, starting with next observation
    virtual real QLearning(real step_size, real gamma)
    {
        //Serror("Not implemented\n");
        return tree.QLearning(step_size, gamma, current_obs, current_reward);
    }

    /// Do Sarsa, starting with next observation
    virtual real Sarsa(real step_size, real gamma, real epsilon)
    {
        Serror("Not implemented\n");
        return tree.Sarsa(step_size, gamma, current_obs, current_reward);
    }

    virtual real ObservationProbability (int act, int x) 
    {
        Serror("Not implemented\n");
        return -1;
    }

    virtual void Reset()
    {
        current_obs = 0;
    }
};

/** A factored predictor RL wrapper for continuous state and action problems. */
template <class T>
class ContinuousTFactoredPredictorRL : public FactoredPredictorRL<Vector, Vector>
{
protected:
    int n_contexts; ///< the number of distinct contexts
    int n_obs; ///< the number of distinct observations
    int n_actions; ///< the number of actions
    T tree; ///< the context tree
    Vector current_obs; ///< the current observation
    real current_reward; ///< the current reward
public:
    ContinuousTFactoredPredictorRL(const Vector& lower_bound_state_action,
								   const Vector& upper_bound_state_action,
								   const Vector& lower_bound_state,
								   const Vector& upper_bound_state,
								   int context_depth,
								   int prediction_depth)
        : n_contexts(lower_bound_state_action.Size()),
          n_obs(lower_bound_state.Size()),
		  n_actions(n_contexts - n_obs),
          tree(2,  context_depth, prediction_depth,
			   lower_bound_state_action,
			   upper_bound_state_action,
			   lower_bound_state,
			   upper_bound_state),
			   current_obs(n_obs)
    {        
		printf ("contexts: %d, obs: %d, actions: %d\n", n_contexts, n_obs, n_actions);
    }


    virtual ~ContinuousTFactoredPredictorRL()
    {
    }
    /* Training and generation */
    /// Observe the (first?) observation.
    virtual real Observe (const Vector& obs) 
    {
        current_obs = obs;
        return 1.0 / (real) n_obs;
    }

    /// Observe current action and next observation and reward.
    ///
    /// As a side-effect, the current observation changes.
    virtual real Observe (const Vector& act, const Vector& obs, real reward) 
    {
		assert(n_obs == obs.Size());
		assert(n_actions == act.Size());

		Vector x(n_actions + n_obs);
		for (int i=0; i<n_obs; ++i) {
			x(i) = current_obs(i);
		}
		for (int i=0; i<n_actions; ++i) {
			x(i + n_obs) = act(i);
		}

		current_obs = obs;
        current_reward = reward;
		
        return tree.Observe(x, obs, reward);
    }
    
    /// Observe current action and next observation
    virtual real QValue (const Vector& act) 
    {
		assert(n_actions == act.Size());

		Vector x(n_actions + n_obs);
		for (int i=0; i<n_obs; ++i) {
			x(i) = current_obs(i);
		}
		for (int i=0; i<n_actions; ++i) {
			x(i + n_obs) = act(i);
		}
        return tree.QValue(x);
    }


    /// Do q-learning, starting with next observation
    virtual real QLearning(real step_size, real gamma)
    {
        //Serror("Not implemented\n");
        return tree.QLearning(step_size, gamma, current_obs, current_reward);
    }

    /// Do q-learning, starting with next observation
    virtual real Sarsa(real step_size, real gamma, real epsilon)
    {
        Serror("Not implemented\n");
        return tree.Sarsa(step_size, gamma, current_obs, current_reward);
    }

    virtual real ObservationProbability (const Vector& act, const Vector& x) 
    {
        Serror("Not implemented\n");
        return -1;
    }

    virtual void Reset()
    {
        current_obs.Clear();
		tree.Reset();
    }
};


/** A factored predictor RL wrapper for continuous state problems */
template <class T>
class ContinuousStateTFactoredPredictorRL : public FactoredPredictorRL<Vector, int>
{
protected:
    int n_contexts; ///< the number of distinct contexts
    int n_obs; ///< the number of distinct observations
    int n_actions; ///< the number of actions
    T tree; ///< the context tree
    Vector current_obs; ///< the current observation, x_{t+1}
	int current_action; ///< the current action, a_t
    real current_reward; ///< the current reward, r_{t+1}
public:
    /**  Construct a predictor.

         n_actions_ the number of actions in the problem
         context_depth: the depth of the conditional context tree
         prediction_depth: the depth of the density esimation tree
         depth_factor: a prior on how the weights should depend on depth
         weight_factor: a prior on the initial value of weights
	*/
    ContinuousStateTFactoredPredictorRL(int n_actions_,
										const Vector& lower_bound_state,
										const Vector& upper_bound_state,
										int context_depth,
                                        int prediction_depth,
                                        real depth_factor,
                                        real weight_factor)
        : n_obs(lower_bound_state.Size()),
		  n_actions(n_actions_),
          tree(n_actions,  context_depth, prediction_depth,
			   lower_bound_state,
			   upper_bound_state,
               depth_factor,
               weight_factor),
		  current_obs(n_obs),
		  current_action(-1)
    {        
		printf (" obs: %d, actions: %d\n", n_obs, n_actions);
    }
	
    /// Destructor
    virtual ~ContinuousStateTFactoredPredictorRL()
    {
    }


    /* --- Training and generation --- */

    /// Observe the (first?) observation.
    virtual real Observe (const Vector& obs) 
    {
        current_obs = obs;
        return 1.0 / (real) n_obs;
    }

    /** Observe current action and next observation and reward.
		
		As a side-effect, the current observation changes.  The model
		assumes that:
		\f[
		x_{t+1} = obs,  r_{t+1} = reward, a_t = a
		\f]
		while the dependencies are Markovian with distributions:
		\f[
		P(x_{t+1} | a_t, x_t), P(r_{t+1} | a_t, x_t).
		\f]

		This has the side-effect of updating the weight values
		of the tree.
	*/
    virtual real Observe (const int& a, const Vector& obs, real reward) 
    {
		assert(n_obs == obs.Size());
		assert(a >= 0 && a < n_actions);

		current_action = a;
		
		real p = tree.Observe(current_obs, current_action, obs, reward);

		current_obs = obs;
        current_reward = reward;
		return p;
    }
    
    /** Get the q-value of actions at current observation.
		
		This must be called after Observe(), of course.
	 */
    virtual real QValue (const int& a) 
    {
		assert(a >= 0 && a < n_actions);
        return tree.QValue(current_obs, a);
    }

    /// Do q-learning, starting with next observation
    virtual real QLearning(real step_size, real gamma)
    {
        //Serror("Not implemented\n");
        return tree.QLearning(step_size, gamma, current_obs, current_reward);
    }

    /// Do q-learning, starting with next observation
    virtual real Sarsa(real step_size, real gamma, real epsilon)
    {
        return tree.Sarsa(epsilon, step_size, gamma, current_obs, current_reward);
    }

    /// Get the probability of a particular observation
    virtual real ObservationProbability (const int& act, const Vector& x) 
    {
        Serror("Not implemented\n");
        return -1;
    }

    /// Reset the model's state, without changing its parameters
    virtual void Reset()
    {
        current_obs.Clear();
		current_action = -1;
		tree.Reset();
    }

    /// Print some statistics about the model to stdout
    virtual void Show()
    {
        tree.Show();
    }
};






#endif
