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

#include "ContextTreeBMC.h"

ContextTreeBMC::Node::Node(int n_branches_,
						int n_outcomes_,
                        std::vector<real>& log_p_)
	: n_branches(n_branches_),
	  n_outcomes(n_outcomes_),
	  depth(0),
	  prev(NULL),
	  next(n_branches),
	  P(n_outcomes), alpha(n_outcomes), prior_alpha(0.5),
      log_p(log_p_)
{
	for (int i=0; i<n_outcomes; ++i) {
		P(i) = 1.0 / (real) n_outcomes;
		alpha(i) = 0;
	}
}

/// Make a node for K symbols at nominal depth d
ContextTreeBMC::Node::Node(ContextTreeBMC::Node* prev_)
	: n_branches(prev_->n_branches),
	  n_outcomes(prev_->n_outcomes),
	  depth(prev_->depth + 1),
	  prev(prev_),
	  next(n_branches),
	  P(n_outcomes),
	  alpha(n_outcomes),
	  prior_alpha(0.5),
      log_p(prev->log_p)
{
	for (int i=0; i<n_branches; ++i) {
		next[i] = NULL;
	}
	for (int i=0; i<n_outcomes; ++i) {
		P(i) = 1.0 / (real) n_outcomes;
		alpha(i) = 0;
	}
    log_p[depth] = 1.0 / (real) n_outcomes;
}

/// make sure to kill all
ContextTreeBMC::Node::~Node()
{
	for (int i=0; i<n_branches; ++i) {
		delete next[i];
	}
}

real ContextTreeBMC::Node::Observe(Ring<int>& history,
								Ring<int>::iterator x,
								int y,
								real probability)
{
    real S = alpha.Sum();
	real Z = 1.0 / (prior_alpha * (real) n_outcomes + S);
	P = (alpha + prior_alpha) * Z;
    log_p[depth] = log(P(y));
	alpha[y]++;

    // P(y | B_k) = P(y | B_k, h_k) P(h_k | B_k) + (1 - P(h_k | B_k)) P(y | B_{k-1})

#if 0
    std::cout << depth << ": P(y|h_k)=" << P[y] 
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif

    // Make sure we have enough observations to justify adding a
    // node. This means at least as many as total outcomes.
    real threshold = (real) n_outcomes; 

    // Go deeper when there has been at least one observations
    // node. 
    //real threshold = 2;

    // Always go deepr, no matter what
    //real threshold = 0; 

    // Go deeper if the context is long enough and the number of
    // observations justifies it.
	if (x != history.end() && S >  threshold) {
		int k = *x;
		++x;
		if (!next[k]) {
			next[k] = new Node(this);
            if (log_p.size() < (uint) depth + 2) {
                log_p.push_back(log (P[y]));
            }
		}
		next[k]->Observe(history, x, y, 0);
	}

	return 0;
}


void ContextTreeBMC::Node::Show()
{
	
    //std::cout << w << " " << depth << "# weight depth\n";
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

int ContextTreeBMC::Node::NChildren()
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

ContextTreeBMC::ContextTreeBMC(int n_branches_, int n_symbols_, int max_depth_)
	: n_branches(n_branches_),
	  n_symbols(n_symbols_),
	  max_depth(max_depth_),
	  history(max_depth),
      w(1),
      log_w(1),
      log_p(1),
      p(1)
{
    log_p[0] = 0;
    log_w[0] = 0;
    p[0] = 0.5;
    w[0] = 1;
	root = new Node(n_branches, n_symbols, log_p);
}

ContextTreeBMC::~ContextTreeBMC()
{
	delete root;
}

real ContextTreeBMC::Observe(int x, int y)
{
    history.push_back(x);
	root->Observe(history, history.begin(), y, 0);

    real log_sum_w = LOG_ZERO;
    //real sum_w = 0;
    for (uint i=0; i<log_w.size(); ++i) {
        log_sum_w = logAdd(log_sum_w, log_w[i]);
        //sum_w += w[i];
    }
    real log_p_x = LOG_ZERO;
    //real p_x = 0;
    for (uint i=0; i<log_w.size(); ++i) {
        log_w[i] -= log_sum_w;
        //w[i] /= sum_w;
        real log_pw = log_p[i] + log_w[i];
        log_p_x = logAdd(log_p_x, log_pw);
        //p_x += w[i] * exp(log_p[i]);
    }

    for (uint i=0; i<log_w.size(); ++i) {
        log_w[i] = log_p[i] + log_w[i] - log_p_x;   
        w[i] = exp(log_w[i]);
        //printf ("%f ", w[i]);
    }
    //    printf ("#WEIGHT\n");

    // Add the new weight and normalise only at the end.
    if (log_p.size() > log_w.size()) {
        uint s = log_w.size();
        log_w.resize(log_p.size());
        w.resize(log_p.size());
        for (uint i=s; i<log_p.size(); ++i) {
            real new_log_w  = - (real) i * log(2);
            real new_w = exp(new_log_w);
            log_w[i] = new_log_w;
            w[i] = new_w;
        }

        real log_sum_w = LOG_ZERO;

        for (uint i=0; i<log_w.size(); ++i) {
            log_sum_w = logAdd(log_sum_w, log_w[i]);
            //sum_w += w[i];
        }
        real sum_w = 0;
        for (uint i=0; i<log_w.size(); ++i) {
            log_w[i] -= log_sum_w;
            //printf("normalised: %f #WEIGHT\n", w[i]);
            w[i] /= sum_w;
        }
    }
    
    return exp(log_p_x);
}

void ContextTreeBMC::Show()
{
    root->Show();
	std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTreeBMC::NChildren()
{
	return root->NChildren();
}
