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

#include "ContextTreeCTW.h"

ContextTreeCTW::Node::Node(int n_branches_,
						int n_outcomes_)
	: n_branches(n_branches_),
	  n_outcomes(n_outcomes_),
	  depth(0),
	  prev(NULL),
	  next(n_branches),
	  P(n_outcomes), alpha(n_outcomes), prior_alpha(0.5),
      w(1), log_w(0), log_w_prior(0)
{
	for (int i=0; i<n_outcomes; ++i) {
		P(i) = 1.0 / (real) n_outcomes;
		alpha(i) = 0;
	}
}

/// Make a node for K symbols at nominal depth d
ContextTreeCTW::Node::Node(ContextTreeCTW::Node* prev_)
	: n_branches(prev_->n_branches),
	  n_outcomes(prev_->n_outcomes),
	  depth(prev_->depth + 1),
	  prev(prev_),
	  next(n_branches),
	  P(n_outcomes),
	  alpha(n_outcomes),
	  prior_alpha(0.5),
	  log_w(0),
	  log_w_prior(prev_->log_w_prior)
{
    w = exp(log_w_prior);
	for (int i=0; i<n_branches; ++i) {
		next[i] = NULL;
	}
	for (int i=0; i<n_outcomes; ++i) {
		P(i) = 1.0 / (real) n_outcomes;
		alpha(i) = 0;
	}

}

/// make sure to kill all
ContextTreeCTW::Node::~Node()
{
	for (int i=0; i<n_branches; ++i) {
		delete next[i];
	}
}

real ContextTreeCTW::Node::Observe(Ring<int>& history,
								Ring<int>::iterator x,
								int y,
								real probability)
{
	real total_probability = 0;
	// calculate probabilities
    real S = alpha.Sum();
	real Z = 1.0 / (prior_alpha * (real) n_outcomes + S);
	P = (alpha + prior_alpha) * Z;
	alpha[y]++;
    total_probability = P[y];

    real threshold = (real) n_outcomes;

	if (x != history.end() && S >  threshold) {
		int k = *x;
		++x;
		if (!next[k]) {
			next[k] = new Node(this);
		}
		total_probability = next[k]->Observe(history, x, y, total_probability);
	}

    total_probability = 0.5 * (P[y] + total_probability);
	return total_probability;
}

void ContextTreeCTW::Node::Show()
{
	
    std::cout << w << " " << depth << "# weight depth\n";
	for (int i=0; i<n_outcomes; ++i) {
		std::cout << alpha[i] << " ";
	}
	for (int k=0; k<n_branches; ++k) {
		if (next[k]) {
			std::cout << "b: " << k << std::endl;
			next[k]->Show();
		}
	}
	std::cout << "<<<<\n";
}

int ContextTreeCTW::Node::NChildren()
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

ContextTreeCTW::ContextTreeCTW(int n_branches_, int n_symbols_, int max_depth_)
	: n_branches(n_branches_),
	  n_symbols(n_symbols_),
	  max_depth(max_depth_),
	  history(max_depth)
{
	root = new Node(n_branches, n_symbols);
}

ContextTreeCTW::~ContextTreeCTW()
{
	delete root;
}

real ContextTreeCTW::Observe(int x, int y)
{
    history.push_back(x);
	return root->Observe(history, history.begin(), y, 0);
}

void ContextTreeCTW::Show()
{
    root->Show();
	std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTreeCTW::NChildren()
{
	return root->NChildren();
}
