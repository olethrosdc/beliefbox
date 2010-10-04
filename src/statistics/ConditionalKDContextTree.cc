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

#include "ConditionalKDContextTree.h"
#include "BetaDistribution.h"
#include <cmath>

#include "Random.h"

ConditionalKDContextTree::Node::Node(ConditionalKDContextTree& tree_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(tree_),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(0),
      prev(NULL),
      next(tree.n_branches),
      w(1),
	  log_w(0),
	  //log_w_prior(0),
      S(0)
{
    assert(lower_bound_x < upper_bound_x);
    //assert(lower_bound_y < upper_bound_y);
	splitting_dimension = ArgMax(upper_bound_x - lower_bound_x);
#ifdef RANDOM_SPLITS
	real phi = 0.1 + 0.8 * urandom();
	mid_point = phi * upper_bound_x[splitting_dimension] + 
		(1 - phi) * lower_bound_x[splitting_dimension];
#else
	mid_point = (upper_bound_x[splitting_dimension] + lower_bound_x[splitting_dimension]) / 2.0;
#endif


	local_density = new ContextTreeKDTree(tree.n_branches, tree.max_depth_cond, tree.lower_bound_y, tree.upper_bound_y);
    int y_dim = tree.upper_bound_y.Size();
	normal_density = new MultivariateNormalUnknownMeanPrecision((tree.upper_bound_y + tree.lower_bound_y)*0.5 , 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
    log_prior_normal = log(0.5);
    prior_normal = 0.5;
}

/// Make a node for K symbols at nominal depth d
ConditionalKDContextTree::Node::Node(ConditionalKDContextTree::Node* prev_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(prev_->tree),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(prev_->depth + 1),
      prev(prev_),
      next(tree.n_branches),
      log_w(-depth * log(2)),
	  //log_w_prior(log(0.9)),
	  //log_w_prior(prev_->log_w_prior - log(2)),
      //log_w_prior(-log(2)),
      S(0)
      //log_w_prior( - log(10))
{
    assert(lower_bound_x < upper_bound_x);
	splitting_dimension = ArgMax(upper_bound_x - lower_bound_x);    
#ifdef RANDOM_SPLITS
	real phi = 0.1 + 0.8 * urandom();
	mid_point = phi * upper_bound_x[splitting_dimension] + 
		(1 - phi) * lower_bound_x[splitting_dimension];
#else
	mid_point = (upper_bound_x[splitting_dimension] + lower_bound_x[splitting_dimension]) / 2.0;
#endif

    w = exp(log_w);
    for (int i=0; i<tree.n_branches; ++i) {
        next[i] = NULL;
    }

	local_density = new ContextTreeKDTree(tree.n_branches,
										  tree.max_depth_cond,
										  tree.lower_bound_y,
										  tree.upper_bound_y);
    int y_dim = tree.upper_bound_y.Size();
	normal_density = new MultivariateNormalUnknownMeanPrecision((tree.upper_bound_y + tree.lower_bound_y)*0.5 , 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
    log_prior_normal = log(0.5);
    prior_normal = 0.5;
}

/// make sure to kill all
ConditionalKDContextTree::Node::~Node()
{
	delete local_density;
    delete normal_density;
    for (int i=0; i<tree.n_branches; ++i) {
        delete next[i];
    }
}

/** Observe new data, adapt parameters.
    
    Each node corresponds to an interval C.
	We only want to model the conditional $\Pr(y \mid x \in C)\f$.
*/
real ConditionalKDContextTree::Node::Observe(Vector& x, Vector& y, real probability)
{
	real total_probability;

    // the local distribution
    //real prior_normal = exp(log_prior_normal);
    real P_tree = local_density->Observe(y);
    real P_normal = normal_density->Observe(y);
    real P_local = prior_normal * P_normal + (1 - prior_normal) * P_tree;
    //log_prior_normal += log(P_normal) - log(P_local);
    prior_normal *= P_normal / P_local;
    //printf ("%f %f = %f -> %f\n", P_tree, P_normal, P_local, prior_normal);
	// Mixture with the previous ones
    //w = exp(log_w_prior + log_w); 
    w = exp(log_w); 
    assert (w >= 0 && w <= 1);
    
	total_probability = P_local * w + (1 - w) * probability;

    // adapt parameters
    //log_w = log(w * P_local / total_probability) - log_w_prior;
    log_w += log(P_local) - log(total_probability);

    // Which interval is the x lying at
    int k;
    if ( x[splitting_dimension] < mid_point) {
        k = 0;
    } else {
        k = 1;
    }
    real threshold = depth; //sqrt((real) depth);
    S++;
#if 0
    std::cout << depth << ": P(y|h_k)=" << P_local
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif
	// Do a forward mixture if there is another node available.
    if ((tree.max_depth==0 || depth < tree.max_depth) && S >  threshold) {
        if (!next[k]) {
            if (k == 0) {
				Vector new_bound_x = upper_bound_x;
				new_bound_x(splitting_dimension) = mid_point;
                next[k] = new Node(this, lower_bound_x, new_bound_x);
            } else {
				Vector new_bound_x = lower_bound_x;
				new_bound_x(splitting_dimension) = mid_point;
                next[k] = new Node(this, new_bound_x, upper_bound_x);
            }
        }
		total_probability = next[k]->Observe(x, y, total_probability);
    }


    return total_probability;
}

/** Recursively obtain \f$P(y | x)\f$.
*/
real ConditionalKDContextTree::Node::pdf(Vector& x, Vector& y, real probability)
{
	real total_probability;

    // the local distribution
    //real P_local = local_density->pdf(y);
    //real prior_normal = exp(log_prior_normal);
    real P_tree = local_density->pdf(y);
    real P_normal = normal_density->pdf(y);
    real P_local = prior_normal * P_normal + (1 - prior_normal) * P_tree;

	// Mix the current one with all previous ones
   ///    w =  exp(log_w_prior + log_w); 
    w =  exp(log_w);
	total_probability = P_local * w + (1 - w) * probability;

#if 0
    std::cout << depth << ": P(y|h_k)=" << P_local
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << ", P(N)= " << prior_normal
              << std::endl;
#endif
    // Which interval is the x lying at
    int k;
    if ( x[splitting_dimension] < mid_point) {
        k = 0;
    } else {
        k = 1;
    }    
	// Do one more mixing step if required
    if (next[k]) {
			total_probability = next[k]->pdf(x, y, total_probability);
    }
    return total_probability;
}


void ConditionalKDContextTree::Node::Show()
{
#if 0
	printf("%d %f #w\n", depth, w);
	for (int k=0; k<tree.n_branches; ++k) {
		if (next[k]) {
			next[k]->Show();
		}
	}
#endif
	return;
}

int ConditionalKDContextTree::Node::NChildren()
{
    int my_children = 0;
    for (int k=0; k<tree.n_branches; ++k) {
        if (next[k]) {
            my_children++;
            my_children += next[k]->NChildren();
        }
    }
    return my_children;
}

ConditionalKDContextTree::ConditionalKDContextTree(int n_branches_,
												   int max_depth_,
												   int max_depth_cond_,
												   Vector& lower_bound_x,
												   Vector& upper_bound_x,
												   Vector& lower_bound_y_,
												   Vector& upper_bound_y_)
    : n_branches(n_branches_),
      max_depth(max_depth_),
      max_depth_cond(max_depth_cond_),
	  lower_bound_y(lower_bound_y_),
	  upper_bound_y(upper_bound_y_)
{
    root = new Node(*this, lower_bound_x, upper_bound_x);
}

ConditionalKDContextTree::~ConditionalKDContextTree()
{
    delete root;
}

/** Obtain \f$\xi_t(y \mid x)\f$ and calculate \f$\xi_{t+1}(w) = \xi_t(w \mid x, y)\f$.
 */
real ConditionalKDContextTree::Observe(Vector& x, Vector& y)
{
    return root->Observe(x, y, 1);
}

/** Obtain \f$\xi_t(y \mid x)\f$ */
real ConditionalKDContextTree::pdf(Vector& x, Vector& y)
{
    return root->pdf(x, y, 1);
}

void ConditionalKDContextTree::Show()
{
    root->Show();
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ConditionalKDContextTree::NChildren()
{
    return root->NChildren();
}
