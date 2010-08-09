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

#ifndef CONTEXT_TREE_RL_H
#define CONTEXT_TREE_RL_H

#include <vector>
#include <list>
#include "real.h"
#include "Vector.h"
#include "Ring.h"
#include "BetaDistribution.h"


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
    
    @see ContextTreeRL
*/
class ContinuousContextTreeRL
{
public:
    // public classes
    struct Node
    {
        int n_branches; ///< number of conditional observations
        int n_outcomes; ///< number of observations to predict
        int depth; ///< depth
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
        Vector P; ///< probability of next symbols
        Vector alpha; ///< parameters of next symbols
        const real prior_alpha; ///< implicit prior value of alpha
        real w; ///< backoff weight
        real log_w; ///< log of w
        real log_w_prior; ///< initial value
        BetaDistribution reward_prior;
        Node(int n_branches_,
             int n_outcomes_);
        Node(Node* prev_);
        ~Node();
        real Q; ///< last Q value of the context
        real w_prod; ///< \f$\prod_k (1 - w_k)\f$
        real context_probability; ///< last probability of the context
        real Observe(Vector& x,
					 int a,
					 Vector& y,
                     real r,
                     real probability,
                     std::list<Node*>& active_contexts);
        real QValue(Vector& x,
					int a,
                    real Q_prev);
        void Show();
        int NChildren();
    };  
    // public methods
    ContinuousContextTreeRL(int n_branches_, 
                  int n_observations,
                  int n_actions,
                  int n_symbols_,
                  int max_depth_= 0);
    ~ContinuousContextTreeRL();
    real Observe(Vector& x, int a, Vector& y, real r);
    void Show();
    int NChildren();
    real QValue(int x);
    real QLearning(real step_size,  real gamma, int observation, real reward);
    real Sarsa(real step_size,  real gamma, int observation, real reward);
protected: 
    int n_branches;
    int n_observations;
    int n_actions;
    int n_symbols;
    int max_depth;
    Node* root;
    Ring<int> history;
    std::list<Node*> active_contexts;
};



#endif
