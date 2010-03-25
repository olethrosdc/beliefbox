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

#include "ContextTreeRealLine.h"
#include "BetaDistribution.h"
#include <cmath>


ContextTreeRealLine::Node::Node(real lower_bound_,
                                real upper_bound_,
                                int n_branches_,
                                int max_depth_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      new_bound((lower_bound + upper_bound)/2),
      n_branches(n_branches_),
      depth(0),
      max_depth(max_depth_),
      prev(NULL),
      next(n_branches),
      alpha(n_branches),
      w(1), log_w(0), log_w_prior(-log(2)),
      w_local(4),
      S(0)
{
    SMART_ASSERT(lower_bound < upper_bound)(lower_bound)(upper_bound);\
    for (int i=0; i<4; ++i) {
        w_local[i] = 0.25;
    }
}

/// Make a node for K symbols at nominal depth d
ContextTreeRealLine::Node::Node(ContextTreeRealLine::Node* prev_,
                                real lower_bound_,
                                real upper_bound_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      n_branches(prev_->n_branches),
      depth(prev_->depth + 1),
      max_depth(prev_->max_depth),
      prev(prev_),
      next(n_branches),
      alpha(n_branches),
      log_w(0),
      log_w_prior(-log(2)),
      w_local(4),
      S(0)
      //log_w_prior( - log(10))
{
    SMART_ASSERT(lower_bound < upper_bound)(lower_bound)(upper_bound)(depth);
    if (isnan(lower_bound) || isnan(upper_bound)) {
        new_bound = 0;
    } else {
        new_bound = (lower_bound + upper_bound) / 2;
    }
    w = exp(log_w_prior);
    for (int i=0; i<4; ++i) {
        w_local[i] = 0.25;
    }
    for (int i=0; i<n_branches; ++i) {
        next[i] = NULL;
    }
    
}

/// make sure to kill all
ContextTreeRealLine::Node::~Node()
{
    for (int i=0; i<n_branches; ++i) {
        delete next[i];
    }
}

/** Observe new data, adapt parameters.
    
    Each node corresponds to an interval C.
    The probability of \f$x\f$ is
    \f[
    P(X | C_k) = U(C_k) w_k + (1 - w_k) [a_{k,0} / n_k P(X | C_{k,0}) + a_{k,1} / n_k P(X | C_{k,1})]
    \f]
    of course, if $X \in C_k^i$ then
    \f[
    P(X | C_k) = U(C_k) w_k + (1 - w_k) [a_{k,i} / n_k P(X | C_{k,i}]
    \f]

    Or in long hand
    \f[
    P(X | C) = P(U | C) P(X | U, C) + [1 - P(U | C)] \sum_{C'} P(X | C') P(C' | C, \not U)
    \f]
    That means that
    \f[
    P(U | X, C) = P(X | U, C) P(U | C) / P(X | C).
    \f]
    
    This algorithm uses the same structure of uniform / non-uniform mixture
    as in Marcus Hutter's work Bayesian Infinite Trees.
    
    However, it is significantly simplified by using an explicit
    online posterior recursion.
*/
real ContextTreeRealLine::Node::Observe(real x,
                                        real probability)
{
    // Which interval is the x lying at
    int k;
    if ( x < new_bound) {
        k = 0;
    } else {
        k = 1;
    }
    
    // the two local distributions
    real P_uniform = 1.0 / (upper_bound - lower_bound);
    Vector P_sub(4);
    P_sub[0] = P_uniform;
    //P_sub[1] = 2 * P_uniform * (1.0 - 2 * fabs((x - new_bound) * P_uniform));
    P_sub[1] = 0;
    P_sub[2]=  2 * P_uniform * (1.0 - fabs((x - lower_bound) * P_uniform));
    P_sub[3] = 2 * P_uniform * (1.0 - fabs((upper_bound - x) * P_uniform));
    Vector w_local_post = P_sub * w_local;
    real P_local = w_local_post.Sum();
    w_local = w_local_post / P_local;


    // probability of recursion
    P =  (1.0 + alpha[k]) / (2.0 + S);


    // adapt parameters
    alpha[k] ++;
    S++;

    real threshold = 2;
    if ((max_depth==0 || depth < max_depth) && S >  threshold) {
        if (!next[k]) {
            if (k == 0) {
                next[k] = new Node(this, lower_bound, new_bound);
            } else {
                next[k] = new Node(this, new_bound, upper_bound);
            }
        }
        P *= next[k]->Observe(x, P);
    } else {
        P *= 2 * P_uniform;
    }

    w = exp(log_w_prior + log_w); 


    //real P_local = w_uni * P_uniform + w_tri * P_tri;
    //w_uni = P_uniform * w_uni / P_local;
    //w_tri = 1 - w_uni;



    //real total_probability = P_uniform * w + (1 - w) * P;
    real total_probability = P_local * w + (1 - w) * P;


    // posterior weight
    //log_w = log(w * P_uniform / total_probability) - log_w_prior;
    log_w = log(w * P_local / total_probability) - log_w_prior;

#if 0
    std::cout << depth << ": P(y|h_k)=" << P
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif

    return total_probability;
}

/** Observe new data, adapt parameters.

    We again ignore the probability given by the previous model,
    since this is a much simpler case.

    We only care about the fact that the previous model posits a
    uniform distribution for the interval of the current model, while
    the current model posits a mixture of two uniform distributions.

*/
real ContextTreeRealLine::Node::pdf(real x,
                                    real probability)
{
    int k;
    if ( x < new_bound) {
        k = 0;
    } else {
        k = 1;
    }
    real P_uniform = 1.0 / (upper_bound - lower_bound);

    P =  (1.0 + alpha[k]) / (2.0 + S);

    real threshold = 1;
    if (S >  threshold && next[k]) {
        P *= next[k]->pdf(x, P);
    } else {
        P *= 2 * P_uniform;
    }

    w = exp(log_w_prior + log_w); 
    //real total_probability = P_uniform * w + (1 - w) * P;
    //real P_local = w_uni * P_uniform + w_tri * P_tri;
    Vector P_sub(4);
    P_sub[0] = P_uniform;
    //P_sub[1] = 2 * P_uniform * (1.0 - 2 * fabs((x - new_bound) * P_uniform));
    //P_sub[0]=  0;
    P_sub[1] = 0;
     P_sub[2]=  2 * P_uniform * (1.0 - fabs((x - lower_bound) * P_uniform));
     P_sub[3] = 2 * P_uniform * (1.0 - fabs((upper_bound - x) * P_uniform));
    Vector w_local_post = P_sub * w_local;
    real P_local = w_local_post.Sum();
    //w_local.print(stdout);    



    real total_probability = P_local * w + (1 - w) * P;

    return total_probability;
}


void ContextTreeRealLine::Node::Show()
{
    std::cout << S << " " << w << " " << depth
              << "[ " << lower_bound << ", " << upper_bound 
              << " ] # obs weight depth\n";
    for (int k=0; k<n_branches; ++k) {
        if (next[k]) {
            std::cout << "b: " << k << std::endl;
            next[k]->Show();
        }
    }
    std::cout << "<<<<\n";
}

int ContextTreeRealLine::Node::NChildren()
{
    int my_children = 0;
    for (int k=0; k<n_branches; ++k) {
        if (next[k]) {
            my_children++;
            my_children += next[k]->NChildren();
        }
    }
    return my_children;
}

ContextTreeRealLine::ContextTreeRealLine(int n_branches_,
                                         int max_depth_,
                                         real lower_bound,
                                         real upper_bound)
    : n_branches(n_branches_),
      max_depth(max_depth_)
{
    root = new Node(lower_bound, upper_bound, n_branches, max_depth);
}

ContextTreeRealLine::~ContextTreeRealLine()
{
    delete root;
}

real ContextTreeRealLine::Observe(real x)
{
    return root->Observe(x, 1);
}


real ContextTreeRealLine::pdf(real x)
{
    return root->pdf(x, 1);
}

void ContextTreeRealLine::Show()
{
    root->Show();
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTreeRealLine::NChildren()
{
    return root->NChildren();
}
