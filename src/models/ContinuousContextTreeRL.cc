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

// NOTE: This should use ConditionalKDContextTree


#include "ContinuousContextTreeRL.h"
#include "BetaDistribution.h"
#include <cmath>

#include "Random.h"

ContinuousContextTreeRL::Node::Node(ContinuousContextTreeRL& tree_,
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
	  log_w_prior(0),
      S(0)
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

	
	local_density = new ContextTreeKDTree(tree.n_branches, tree.max_depth_cond, tree.lower_bound_y, tree.upper_bound_y);
    int y_dim = tree.upper_bound_y.Size();
	Vector mean  = (tree.upper_bound_y + tree.lower_bound_y)*0.5;
	printf("y_dim: %d %d- ", y_dim, tree.lower_bound_y.Size()); 
	mean.print(stdout);
	tree.lower_bound_y.print(stdout);
	tree.upper_bound_y.print(stdout);
	normal_density = new MultivariateNormalUnknownMeanPrecision(mean, 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
    log_prior_normal = log(0.5);
    prior_normal = 0.5;
}

/// Make a node for K symbols at nominal depth d
ContinuousContextTreeRL::Node::Node(ContinuousContextTreeRL::Node* prev_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(prev_->tree),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(prev_->depth + 1),
      prev(prev_),
      next(tree.n_branches),
      log_w(0),
	  //log_w_prior(log(0.9)),
	  log_w_prior(prev_->log_w_prior - log(2)),
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

    w = exp(log_w_prior);
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
ContinuousContextTreeRL::Node::~Node()
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

	x is the curent context
	y is the next observation
	r is the next reward
*/
real ContinuousContextTreeRL::Node::Observe(Vector& x, Vector& y, real reward, real probability, ContextList& active_contexts)
{
	active_contexts.push_back(this);
	real total_probability;

    // the local distribution
    //real prior_normal = exp(log_prior_normal);
    real P_tree = local_density->Observe(y);
    real P_normal = normal_density->Observe(y);
    real P_local = prior_normal * P_normal + (1 - prior_normal) * P_tree;
    //log_prior_normal += log(P_normal) - log(P_local);
    prior_normal *= P_normal / P_local;
    //printf ("%f %f = %f -> %f\n", P_tree, P_normal, P_local, prior_normal);

	real p_reward = reward_prior.Observe(reward);
	P_local *= p_reward;
 	// Mixture with the previous ones
    w = exp(log_w_prior + log_w); 
	total_probability = P_local * w + (1 - w) * probability;

    // adapt parameters
    log_w = log(w * P_local / total_probability) - log_w_prior;

	w_prod = 1; ///< auxilliary calculation

    // Which interval is the x lying at
    int k;
    if ( x[splitting_dimension] < mid_point) {
        k = 0;
    } else {
        k = 1;
    }

    real threshold = sqrt((real) depth);
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
		total_probability = next[k]->Observe(x, y, reward, total_probability, active_contexts);
		w_prod = next[k]->w_prod; 
        assert(!isnan(total_probability));
        assert(!isnan(w_prod));
    }


	// Auxilliary calculation for context
    context_probability = w * w_prod;
    w_prod *= (1 - w);

    assert(!isnan(w_prod));
    assert(!isnan(context_probability));

    return total_probability;
}

real ContinuousContextTreeRL::Node::QValue(Vector& x,
										   real Q_prev)
{
	
    w = exp(log_w_prior + log_w); 
    real Q_next = Q * w + (1 - w) * Q_prev;
	printf ("%f -(%f,%f)-> %f\n", Q_prev, Q, w, Q_next);
    if (isnan(Q_prev)) {
        fprintf(stderr, "Warning: at depth %d, Q_prev is nan\n", depth);
        Q_prev = 0;
    }
    if (isnan(Q)) {
        fprintf(stderr, "Warning: at depth %d, Q is nan\n", depth);
        Q = 0;
    }
    if (isnan(Q_next)) {
        fprintf(stderr, "Warning: at depth %d, Q_next is nan\n", depth);
        Q_next = 0;
    }
    if (isnan(w)) {
        fprintf(stderr, "Warning: at depth %d, w is nan\n", depth);
        w = 0.5;
    }

    int k;
    if ( x[splitting_dimension] < mid_point) {
        k = 0;
    } else {
        k = 1;
    }
	
	if (next[k]) {
        Q_next = next[k]->QValue(x, Q_next);
    } 

    return Q_next;
}

/** Perform one Q-learning step.

	The step-size is the amount by which to change the Q-values.
	gamma is the discount factor.
	y is the next observation
	r is the next reward.

	The function returns the TD-error.
 */
real ContinuousContextTreeRL::QLearning(real step_size, real gamma, Vector& y, real reward)
{
	real Q_prev = root->QValue(current_state_action, 0);

    if (isnan(Q_prev)) {
        Q_prev = 0;
        fprintf(stderr, "Warning: Q_prev is nan\n");
    }

    real max_Q = -INF;    
	Vector x(n_states + n_actions);
	for (int i=0; i<n_states; ++i) {
		x(i) = y(i);
	}
    for (int a = 0; a < n_actions; ++a) {
		for (int i=0; i<n_actions; ++i) {
			if (a==i) {
				x(i) = 1;
			} else {
				x(i) = 0;
			}
		}
        real Q_x = QValue(x);
        if (Q_x > max_Q) {
            max_Q = Q_x;
        }

		printf("Q[%d] = %f\n", a, Q_x);
    }

    real td_err = 0;
    real dQ_i = reward + gamma * max_Q - Q_prev; 
    real p = 0;
    for (std::list<Node*>::iterator i = active_contexts.begin();
         i != active_contexts.end();
         ++i) {
        real p_i = (*i)->context_probability;
        p += p_i;
        real delta = p_i * dQ_i;  // proper way to do it..
        //real delta = dQ_i; 
        //real delta = reward + gamma * max_Q - (*i)->Q; // Alternative approach.
        (*i)->Q += step_size * delta;
        //printf ("%f * %f = %f ->  %f\n", p_i, dQ_i, delta, (*i)->Q);
        //td_err += fabs(delta);
    }
    td_err = fabs(dQ_i);
    printf ("# max_Q:%f Q:%f r:%f TD:%f, (%f %f %f)\n",
			max_Q, Q_prev, reward, td_err, step_size, dQ_i, p);
    return td_err;
}




void ContinuousContextTreeRL::Node::Show()
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

int ContinuousContextTreeRL::Node::NChildren()
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

ContinuousContextTreeRL::ContinuousContextTreeRL(int n_branches_,
												   int max_depth_,
												   int max_depth_cond_,
												   Vector& lower_bound_x,
												 Vector& upper_bound_x,
												   Vector& lower_bound_y_,
												   Vector& upper_bound_y_)
    : 	n_branches(n_branches_),
		max_depth(max_depth_),
		max_depth_cond(max_depth_cond_),
		lower_bound_y(lower_bound_y_),
		upper_bound_y(upper_bound_y_)
{
    root = new Node(*this, lower_bound_x, upper_bound_x);
	n_actions= lower_bound_y.Size();
	n_states = lower_bound_x.Size() - n_actions;
	assert(n_actions > 0);
	assert(n_states > 0);
	current_state_action.Resize(n_actions + n_states);
}

ContinuousContextTreeRL::~ContinuousContextTreeRL()
{
    delete root;
}

/** Obtain \f$\xi_t(y \mid x)\f$ and calculate \f$\xi_{t+1}(w) = \xi_t(w \mid x, y)\f$.
 */
real ContinuousContextTreeRL::Observe(Vector& x, Vector& y, real r)
{
	active_contexts.clear();
    return root->Observe(x, y, r, 1, active_contexts);
}

/** Obtain \f$\xi_t(y \mid x)\f$ and calculate \f$\xi_{t+1}(w) = \xi_t(w \mid x, y)\f$.
 */
real ContinuousContextTreeRL::QValue(Vector& x)
{
    return root->QValue(x, -9999);
}

void ContinuousContextTreeRL::Show()
{
    root->Show();
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContinuousContextTreeRL::NChildren()
{
    return root->NChildren();
}

