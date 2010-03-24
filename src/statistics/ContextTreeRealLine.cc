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
      w(1), log_w(0), log_w_prior(0),
      S(0)
{
    SMART_ASSERT(lower_bound < upper_bound)(lower_bound)(upper_bound);
}

/// Make a node for K symbols at nominal depth d
ContextTreeRealLine::Node::Node(ContextTreeRealLine::Node* prev_,
                                real lower_bound_,
                                real upper_bound_)
	: lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      new_bound((lower_bound + upper_bound)/2),
      n_branches(prev_->n_branches),
	  depth(prev_->depth + 1),
      max_depth(prev_->max_depth),
	  prev(prev_),
	  next(n_branches),
      alpha(n_branches),
	  log_w(0),
	  log_w_prior(prev_->log_w_prior - log(2)),
      S(0)
	  //log_w_prior( - log(10))
{
    SMART_ASSERT(lower_bound < upper_bound)(lower_bound)(upper_bound);
    w = exp(log_w_prior);
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

    We basically ignore the probability given by the previous model,
    since this is a much simpler case.

    We only care about the fact that the previous model posits a
    uniform distribution for the interval of the current model, while
    the current model posits a mixture of two uniform distributions.

    Thus, the only thing that matters (for estimating the weights) is
    the ratio of uniform to non-uniform probabilities.

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

    real Z = 1.0 / (upper_bound - lower_bound);
    real P_uniform =  Z * 0.5;
    P = Z * (1.0 + alpha[k]) / (2.0 + S);

    alpha[k] ++;

    w = exp(log_w_prior + log_w); 

    real total_probability = P * w + (1 - w) * P_uniform;

#if 0
    std::cout << depth << ": P(y|h_k)=" << P
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif
    log_w = log(w * P /total_probability) - log_w_prior;

    // Go deeper when there has been at least one observations
    // at the node. 
    real threshold = 1;
	if ((max_depth==0 || depth < max_depth) && S >  threshold) {
		if (!next[k]) {
            if (k == 0) {
                next[k] = new Node(this, lower_bound, new_bound);
            } else {
                next[k] = new Node(this, new_bound, upper_bound);
            }
		}
		total_probability = next[k]->Observe(x, total_probability);
	}

    S++;
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
    
    real Z = 1.0 / (upper_bound - lower_bound);
    real P_uniform =  Z * 0.5;
    P = Z * (1.0 + alpha[k]) / (2.0 + S);

    w = exp(log_w_prior + log_w); 
    real total_probability =  P * w + (1 - w) * P_uniform;
    // Go deeper when there has been at least one observations
    // at the node. 
    real threshold = 1;
	if (S >  threshold) {
		if (next[k]) {
            total_probability = next[k]->pdf(x, total_probability);
        }
	}
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

ContextTreeRealLine::ContextTreeRealLine(int n_branches_, int max_depth_)
	: n_branches(n_branches_),
	  max_depth(max_depth_)
{
	root = new Node(0, 1, n_branches, max_depth);
}

ContextTreeRealLine::~ContextTreeRealLine()
{
	delete root;
}

real ContextTreeRealLine::Observe(real x)
{
	return root->Observe(x, 0);
}


real ContextTreeRealLine::pdf(real x)
{
	return root->pdf(x, 0);
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
