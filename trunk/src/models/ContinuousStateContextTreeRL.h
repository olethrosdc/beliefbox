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

#ifndef CONTINUOUS_STATE_CONTEXT_TREE_RL_H
#define CONTINUOUS_STATE_CONTEXT_TREE_RL_H

#include <vector>
#include <list>
#include "real.h"
#include "Vector.h"
#include "BetaDistribution.h"
#include "NormalDistribution.h"
#include "ContextTreeKDTree.h"

/** A continuous context tree implementation.
    
    This tree partitions an \f$nf\$-dimensional space.  It is very
    similar to ContextTreeRL, the only difference being that the
    partition is done in the state-space.

	There is a different tree for each action.
*/
class ContinuousStateContextTreeRL
{
public:
	struct Node;
    typedef std::list<Node*> ContextList;
    // public classes
    struct Node
    {
		ContinuousStateContextTreeRL& tree; ///< the tree
        Vector lower_bound_x; ///< looks at x > lower_bound
        Vector upper_bound_x; ///< looks at x < upper_bound
        real mid_point; ///< how to split
		int splitting_dimension; ///< dimension on which to do the split.
        const int depth; ///< depth of the node
        real log_prior_normal;
        real prior_normal;
		ContextTreeKDTree* local_density; ///< local density estimator
		MultivariateNormalUnknownMeanPrecision* normal_density; ///< local density estimator
		NormalUnknownMeanPrecision reward_prior; ///< local density estimator
        real mean_reward; ///< mean reward
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
        real log_w; ///< log of w
        real w; ///< backoff weight
		real Q; ///< Q-value
		real w_prod; ///< \f$\prod_k (1 - w_k)\f$
        real context_probability; ///< last probability of the context

        Node(ContinuousStateContextTreeRL& tree_,
			 Vector& lower_bound_x_,
			 Vector& upper_bound_x_);
        Node(Node* prev_, 
			 Vector& lower_bound_x, Vector& upper_bound_x);
        ~Node();
        real Observe(Vector& x, Vector& y, real reward, real probability, ContextList& active_contexts);
		real QValue(Vector& state, real Q_prev);
        virtual void Show();
        int NChildren();    
        int S;
    };

    // public methods
    ContinuousStateContextTreeRL(int n_actions_,
								 int max_depth_,
								 int max_depth_cond_,
								 Vector& lower_bound_x, Vector& upper_bound_x,
                                 real depth_factor_,
                                 real weight_factor_);
    ~ContinuousStateContextTreeRL();
    real Observe(Vector& x, int a, Vector& y, real r);
    //real pdf(Vector& x, Vector& y);
	real QValue(Vector& x);
	real QValue(Vector& x, int a);
	real QLearning(real step_size, real gamma, Vector& y, real reward);
	real ValueIteration();
	real Sarsa(real epsilon, real step_size, real gamma, Vector& y, real reward);
	void Reset()
	{
		current_action = -1;
		active_contexts.clear();
	}
    virtual void Show();
    int NChildren();
protected: 
	int n_states;
	int n_actions;
    int n_branches;
    int max_depth;
    int max_depth_cond;
	Vector lower_bound_x;
	Vector upper_bound_x;
    real depth_factor;
    real weight_factor;
	std::vector<Node*> root;
	Vector current_state; ///< current state 
	int current_action; ///< current action
	std::list<Node*> active_contexts;
};

#endif
