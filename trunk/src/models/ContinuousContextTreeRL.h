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

#ifndef CONTINUOUS_CONTEXT_TREE_RL_H
#define CONTINUOUS_CONTEXT_TREE_RL_H

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

    The simplest models are:

    1. Each context model predicts only observations within that level.
    So, \f$\Pr(z_{t+1} \in c' \mid z_t \in c)\f$, with \f$c, c' \in C_k\f$.
    
    2. Each context model predicts next contexts in general
    So, \f$\Pr(z_{t+1} \in c' \mid z_t \in c)\f$, w'ith \f$c, c' \in C\f$.

    3. There is a full density estimation.

	The observations are the concatenations of actions and observations.

    @see ContextTreeRL, ConditionalKDContextTree
*/
class ContinuousContextTreeRL
{
public:
	struct Node;
    typedef std::list<Node*> ContextList;
    // public classes
    struct Node
    {
		ContinuousContextTreeRL& tree; ///< the tree
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
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
        real w; ///< backoff weight
        real log_w; ///< log of w
        real log_w_prior; ///< initial value
		real Q; ///< Q-value
		real w_prod; ///< \f$\prod_k (1 - w_k)\f$
        real context_probability; ///< last probability of the context

        Node(ContinuousContextTreeRL& tree_,
			 Vector& lower_bound_x_,
			 Vector& upper_bound_x_);
        Node(Node* prev_, 
			 Vector& lower_bound_x, Vector& upper_bound_x);
        ~Node();
        real Observe(Vector& x, Vector& y, real reward, real probability, ContextList& active_contexts);
		real QValue(Vector& state_action, real Q_prev);
        void Show();
        int NChildren();    
        int S;
    };

    // public methods
    ContinuousContextTreeRL(int n_branches_,
							int max_depth_,
                             int max_depth_cond_,
							 Vector& lower_bound_x, Vector& upper_bound_x,
							 Vector& lower_bound_y, Vector& upper_bound_y);
    ~ContinuousContextTreeRL();
    real Observe(Vector& x, Vector& y, real r);
    //real pdf(Vector& x, Vector& y);
	real QValue(Vector& x);
	real QLearning(real step_size, real gamma, Vector& y, real reward);
	real Sarsa(real epsilon, real step_size, real gamma, Vector& y, real reward) {
        Serror("Not implemented\n");
        return 0;
    }
	void Reset()
	{
		active_contexts.clear();
	}
    void Show();
    int NChildren();
protected: 
	int n_states;
	int n_actions;
    int n_branches;
    int max_depth;
    int max_depth_cond;
	Vector lower_bound_y;
	Vector upper_bound_y;
    Node* root;
	Vector current_state_action; ///< current state and action pair
	std::list<Node*> active_contexts;
};

#endif
