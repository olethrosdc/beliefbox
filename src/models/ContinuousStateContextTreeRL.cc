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


#include "ContinuousStateContextTreeRL.h"
#include "BetaDistribution.h"
#include <cmath>

#include "Random.h"
#define INITIAL_Q_VALUE 0.0

ContinuousStateContextTreeRL::Node::Node(ContinuousStateContextTreeRL& tree_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(tree_),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(0),
      mean_reward(0),
      prev(NULL),
      next(tree.n_branches),
      w(1),
	  log_w(0),
	  log_w_prior(0),
	  Q(INITIAL_Q_VALUE),
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

	
	local_density = new ContextTreeKDTree(tree.n_branches, tree.max_depth_cond, tree.lower_bound_x, tree.upper_bound_x);
    int y_dim = tree.upper_bound_x.Size();
	Vector mean  = (tree.upper_bound_x + tree.lower_bound_x)*0.5;
	printf("y_dim: %d %d- ", y_dim, tree.lower_bound_x.Size()); 
	mean.print(stdout);
	tree.lower_bound_x.print(stdout);
	tree.upper_bound_x.print(stdout);
	normal_density = new MultivariateNormalUnknownMeanPrecision(mean, 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
    log_prior_normal = log(0.5);
    prior_normal = 0.5;
}

/// Make a node for K symbols at nominal depth d
ContinuousStateContextTreeRL::Node::Node(ContinuousStateContextTreeRL::Node* prev_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(prev_->tree),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(prev_->depth + 1),
      mean_reward(0),
      prev(prev_),
      next(tree.n_branches),
      log_w(0),
	  //log_w_prior(log(0.9)),
	  log_w_prior(prev_->log_w_prior - log(2)),
      Q(INITIAL_Q_VALUE),
      //Q(prev_->Q),
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
										  tree.lower_bound_x,
										  tree.upper_bound_x);
    int y_dim = tree.upper_bound_x.Size();
	normal_density = new MultivariateNormalUnknownMeanPrecision((tree.upper_bound_x + tree.lower_bound_x)*0.5 , 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
    log_prior_normal = log(0.5);
    prior_normal = 0.5;
}

/// make sure to kill all
ContinuousStateContextTreeRL::Node::~Node()
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
real ContinuousStateContextTreeRL::Node::Observe(Vector& x, Vector& y, real reward, real probability, ContextList& active_contexts)
{
	active_contexts.push_back(this);
	real total_probability;
    real fudge = MIN_PRECISION;
    // the local distribution
    //real prior_normal = exp(log_prior_normal);
    real P_tree = local_density->Observe(y);
    real P_normal = normal_density->Observe(y);
    real P_local =  prior_normal * P_normal + (1 - prior_normal) * P_tree;
    //log_prior_normal += log(P_normal) - log(P_local);
    if (P_local > fudge) {
        prior_normal *= P_normal / P_local;
    } else {
        P_local = fudge;
    }
    //printf ("%f %f = %f -> %f\n", P_tree, P_normal, P_local, prior_normal);

	real p_reward = reward_prior.Observe(reward);
    real mean_reward = reward_prior.getMean();
	P_local *= p_reward;

    if (P_local < fudge) {
        P_local = fudge;
    }
 	// Mixture with the previous ones
    w = exp(log_w_prior + log_w); 
    if (w > 1) {
        w = 1;
    } else if (w < 0) {
        w = 0;
    }
	total_probability = P_local * w + (1 - w) * probability;

    // adapt parameters
    //log_w = log(w * P_local / total_probability) - log_w_prior;
    log_w =  log_w + log(P_local) - log(total_probability);
    assert(!isnan(log_w));

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

real ContinuousStateContextTreeRL::Node::QValue(Vector& x,
										   real Q_prev)
{
	
    w = exp(log_w_prior + log_w); 
    real Q_next = Q * w + (1 - w) * Q_prev;
	//printf ("%f -(%f,%f)-> %f\n", Q_prev, Q, w, Q_next);
    assert (!isnan(Q_prev));
    assert (!isnan(Q));
    assert (!isnan(Q_next));
    assert (!isnan(w));

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

	We wish to move Q(s,a) towards r + \max_a' Q(s',a').

	The function returns the TD-error.
 */
real ContinuousStateContextTreeRL::QLearning(real step_size, real gamma, Vector& y, real reward)
{
	if (current_action == -1) {
		return 0;
	}
    real max_Q = QValue(y);

	real Q_prev = QValue(current_state, current_action);

    assert (!isnan(Q_prev));

    real td_err = 0;
    real dQ_i = reward + gamma * max_Q - Q_prev; 
    real p = 0;
    for (std::list<Node*>::iterator i = active_contexts.begin();
         i != active_contexts.end();
         ++i) {
        real p_i = (*i)->context_probability;
        p += p_i;
        real delta = p_i * dQ_i; 
        (*i)->Q += step_size * delta;
    }
    td_err = fabs(dQ_i);
    //printf ("# max_Q:%f Q:%f r:%f TD:%f, (%f %f %f)\n",
	//max_Q, Q_prev, reward, td_err, step_size, dQ_i, p);
    return td_err;
}

/** Perform value iteration.
    
    Randomly select an x.
    
    Select a context that matches.
    
    Update value of said context.

    NOT IMPLEMENTED
 */
real ContinuousStateContextTreeRL::ValueIteration()
{
    

    return 0;
}


void ContinuousStateContextTreeRL::Node::Show()
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

int ContinuousStateContextTreeRL::Node::NChildren()
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

ContinuousStateContextTreeRL::ContinuousStateContextTreeRL(int n_actions_,
														   int max_depth_,
														   int max_depth_cond_,
														   Vector& lower_bound_x_,
														   Vector& upper_bound_x_)
    : 	n_states(lower_bound_x_.Size()),
		n_actions(n_actions_),
		n_branches(2),
		max_depth(max_depth_),
		max_depth_cond(max_depth_cond_),
		lower_bound_x(lower_bound_x_),
		upper_bound_x(upper_bound_x_),
		root(n_actions)
{
    for (int i=0; i<n_actions; ++i) {
		root[i] = new Node(*this, lower_bound_x, upper_bound_x);
	}
	assert(n_states == upper_bound_x.Size());
	assert(n_actions > 0);
	assert(n_states > 0);
	assert(n_branches > 0);
	current_state.Resize(n_states);
	current_action = -1;
}

ContinuousStateContextTreeRL::~ContinuousStateContextTreeRL()
{
	for (int i=0; i<n_actions; ++i) {
		delete root[i];
	}
}

/** Observe a transition and update the model and active contexts.

	x is the current state
	a is the current action
	y is the next state
	r is the next reward
 */
real ContinuousStateContextTreeRL::Observe(Vector& x, int a, Vector& y, real r)
{
	assert(a >= 0 && a < n_actions);
	active_contexts.clear();
	current_state = x;
	current_action = a;
    return root[a]->Observe(x, y, r, 1, active_contexts);
}

/** Calculate \f$Q(x,a)\f$
 */
real ContinuousStateContextTreeRL::QValue(Vector& x, int a)
{
    return root[a]->QValue(x, -9999);
}

/** Calculate \f$\max_a Q(x,a)\f$
 */
real ContinuousStateContextTreeRL::QValue(Vector& x)
{
	real Q_max = -INF;
	int a_max = -1;
	for (int i=0; i<n_actions; ++i) { 
		real Q_i = root[i]->QValue(x, -9999);
		if (Q_i > Q_max) {
			Q_max = Q_i;
			a_max = i;
		}
	}

    return Q_max;
}

/// Show some statistics about the tree
void ContinuousStateContextTreeRL::Show()
{
	for (int i=0; i<n_actions; ++i) {
		root[i]->Show();
	}
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

/// Get the number of children.
int ContinuousStateContextTreeRL::NChildren()
{
	int n_children = 0;
	for (int i=0; i<n_actions; ++i) {
		n_children += root[i]->NChildren();
	}
    return n_children;
}

