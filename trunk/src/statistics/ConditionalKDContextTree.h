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

#ifndef CONDITIONAL_KD_CONTEXT_TREE_H
#define CONDITIONAL_KD_CONTEXT_TREE_H

#include <vector>
#include "real.h"
#include "Vector.h"
#include "Ring.h"
#include "ContextTreeKDTree.h"

/** Context tree non-parametric conditional density estimation on \f$R^n \times R^m\f$.

    This is a generalisation of the binary tree on [a,b] implemented
	in ContextTreeRealLine to a KD-tree. For any sample \f$x^t \sim
	D^t\f$, with $x_i \in R^n\f$ it estimates the measure
	\f[
	P^t(w) = P(w \mid x^t) \propto P_w(x^t)P^0(w)
	\f]
	and so the marginal
	\f[
	P(x_{t+1} \mid x^t) = \sum_w P_w(x_t) P(w \mid x^t).
	\f]
	The distribution \f$P^t(w)\f$ is simply a product distribution.

	The model can also be used for estimating conditional densities, directly.
	However, this is perhaps not a good idea. 
 */
class ConditionalKDContextTree
{
public:
    // public classes
    struct Node
    {
		ConditionalKDContextTree& tree; ///< the tree
        Vector lower_bound_x; ///< looks at x > lower_bound
        Vector upper_bound_x; ///< looks at x < upper_bound
        real mid_point; ///< how to split
		int splitting_dimension; ///< dimension on which to do the split.
        const int depth; ///< depth of the node
		ContextTreeKDTree* local_density; ///< local density estimator
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
        real w; ///< backoff weight
        real log_w; ///< log of w
        real log_w_prior; ///< initial value

        Node(ConditionalKDContextTree& tree_,
			 Vector& lower_bound_x_,
			 Vector& upper_bound_x_);
        Node(Node* prev_, 
			 Vector& lower_bound_x, Vector& upper_bound_x);
        ~Node();
        real Observe(Vector& x, Vector& y, real probability);
        real pdf(Vector& x, Vector& y, real probability);
        void Show();
        int NChildren();    
        int S;
    };
    
    // public methods
    ConditionalKDContextTree(int n_branches_, int max_depth_, int max_depth_cond_,
							 Vector& lower_bound_x, Vector& upper_bound_x,
							 Vector& lower_bound_y, Vector& upper_bound_y);
    ~ConditionalKDContextTree();
    real Observe(Vector& x, Vector& y);
    real pdf(Vector& x, Vector& y);
    void Show();
    int NChildren();
protected: 
    int n_branches;
    int max_depth;
    int max_depth_cond;
	Vector lower_bound_y;
	Vector upper_bound_y;
    Node* root;
};



#endif
