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

#include "ContextTreePPM.h"

ContextTreePPM::Node::Node(int n_branches_,
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
ContextTreePPM::Node::Node(ContextTreePPM::Node* prev_)
	: n_branches(prev_->n_branches),
	  n_outcomes(prev_->n_outcomes),
	  depth(prev_->depth + 1),
	  prev(prev_),
	  next(n_branches),
	  P(n_outcomes),
	  alpha(n_outcomes),
	  prior_alpha(0.5),
	  log_w(0),
	  log_w_prior(prev_->log_w_prior - log(2))
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
ContextTreePPM::Node::~Node()
{
	for (int i=0; i<n_branches; ++i) {
		delete next[i];
	}
}

/** Update nodes using PPM-C.

    In this variant, the probability of a symbol is
    \f[
    P(x|c) = \alpha(x;c) / [Z(c) + \sum_{y \in X} \alpha(y;c)],
    \f]
    where \f$\alpha(x;c)\f$ is the number of times symbol \f$x\f$
    appeared in context \f$c\f$ and \f$Z\f$ is the number of non-zero
    components of \f$\alpha\f$. The probability of escape is
    \f[
    P(s | c) = Z(c) / [Z(c) + \sum_{y \in X} \alpha(y;c)].
    \f]
    The total probability is recursively expressed as
    \f[
    P(x | x^t) = P(x | c_k) I(x | c_k) + [1 - I(x | c_k)] P(s | c_{k-1}),
    \f]
    where \f$I(x|c_k) = 1 \f$ if  \f$\alpha(x;c_k) > 0 \f$ and 0 otherwise.

    To implement this, we must use a forward mechanism. The algorithm actually
    needs a hack because sometimes the escape probability is 0.
 */
real ContextTreePPM::Node::Observe(Ring<int>& history,
								Ring<int>::iterator x,
								int y,
                                   real probability,
                                   Vector* P_prev)
{
	real total_probability = 0;

	// first calculate P_k(x)
    real S = alpha.Sum();
    real Z = 0;
    for (int i=0; i<n_outcomes; ++i) {
        if (alpha(i) > 0) {
            Z++;
        }
    }
    real iSZ;
    if (S > 0) {
#if 0
        // normal mechanism
        iSZ = 1.0 / (S + Z);
        P = alpha * iSZ;
#else
        // exclusion
        iSZ = 1.0 / (S + Z - 1);
        P = alpha * iSZ;
        //P /= P.Sum();
#endif
    } else {
        iSZ = 1.0 / (real) n_outcomes;
        P = (alpha + 1) * iSZ;
        //P /= P.Sum();
    }

    // now calculate P(s | x)
    real escape = Z * iSZ;
    
    real p_uniform = 1.0  / (real) n_outcomes;
    for (int i=0; i<n_outcomes; ++i) {

        if (alpha[i] == 0) {
            real p = p_uniform;
            if (P_prev) {
                p = (*P_prev)(i);
            }
            if (depth) {
                if (escape)  {
                    P[i] = escape * p;
                }  else {
                    P[i] = p;
                }
            } else {
                P[i] = p_uniform;
            }
        }
    }
    P /= P.Sum();
    total_probability = P[y];

    //printf("P = %f, sum P = %f\n", P[y], P.Sum());

    //printf("%f -> %f (%f) -> %f\n", 
    //       probability, P[y], escape, total_probability);
    // update model
	alpha[y]++;

    real threshold = 1;//(real) n_outcomes;
    real S_next = 0;
	if (x != history.end() && S >  threshold) {
		int k = *x;
		++x;
		if (!next[k]) {
			next[k] = new Node(this);
		} else {
            S_next = next[k]->TotalObservations();
        }
		total_probability = next[k]->Observe(history, x, y, total_probability, &P);
	}
    //real ratio = (1 + S) / (1 + S + S_next);
    //total_probability = (1 - ratio) * P[y] + (ratio) * total_probability;
	return total_probability;
}


void ContextTreePPM::Node::Show()
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

int ContextTreePPM::Node::NChildren()
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

ContextTreePPM::ContextTreePPM(int n_branches_, int n_symbols_, int max_depth_)
	: n_branches(n_branches_),
	  n_symbols(n_symbols_),
	  max_depth(max_depth_),
	  history(max_depth)
{
	root = new Node(n_branches, n_symbols);
}

ContextTreePPM::~ContextTreePPM()
{
	delete root;
}

real ContextTreePPM::Observe(int x, int y)
{
    history.push_back(x);
	return root->Observe(history, history.begin(), y, 0, NULL);
}

void ContextTreePPM::Show()
{
    root->Show();
	std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTreePPM::NChildren()
{
	return root->NChildren();
}
