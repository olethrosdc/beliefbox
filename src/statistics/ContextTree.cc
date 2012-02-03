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

#include "ContextTree.h"

ContextTree::Node::Node(int n_branches_,
						int n_outcomes_)
	: N_obs(0),
      n_branches(n_branches_),
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
ContextTree::Node::Node(ContextTree::Node* prev_)
	: N_obs(0),
      n_branches(prev_->n_branches),
	  n_outcomes(prev_->n_outcomes),
	  depth(prev_->depth + 1),
	  prev(prev_),
	  next(n_branches),
	  P(n_outcomes),
	  alpha(n_outcomes),
	  prior_alpha(0.5),
	  log_w(0),
	  log_w_prior(prev_->log_w_prior - log(2))
	  //log_w_prior( - log(10))
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
ContextTree::Node::~Node()
{
	for (int i=0; i<n_branches; ++i) {
		delete next[i];
	}
}

real ContextTree::Node::Observe(Ring<int>& history,
								Ring<int>::iterator x,
								int y,
								real probability)
{
	real total_probability = 0;
	// calculate probabilities

    // Standard
#if 0
    real S = (real) N_obs; //xalpha.Sum();
    real Z = 1.0 / (prior_alpha * (real) n_outcomes + S);
    P = (alpha + prior_alpha) * Z;
    
    //printf("%f\n", P.Sum());//P.print(stdout);
    //P /= P.Sum();
    //P[y] = (alpha[y] + prior_alpha) * Z;
#endif


#if 1
    // P-BVMM
    // use a modified PPM -- best for many outcomes
    real S = alpha.Sum(); // total number of observations
    real Z = 0; // total prior mass for observed values
    int zeros = 0; // number of unobserved values
    for (int i=0; i<n_outcomes; ++i) {
        if (alpha(i) > 0) {
            Z += prior_alpha; 
        } else {
            zeros++;
        }
    }
    if (zeros == n_outcomes) {
        // case 1: observed nothing yet
        for (int i=0; i<n_outcomes; ++i) {
            if (alpha(i) == 0) {
                P(i) = 1.0 / (real) n_outcomes;
            }
        }
    } else {
        // case 2: observed some things
        real iSZ = 1.0 / (S + Z); // S + Z is the total alpha for observed values
        P = (alpha) * iSZ; // Pr of seeing the i-th symbol if you know you are going to see a previsouly seen symbol
        
        real P_rest = 1.0 - P.Sum(); // Remaining mass
        
        if (zeros > 0) {
            for (int i=0; i<n_outcomes; ++i) {
                if (alpha(i) == 0) {
                    P(i) = P_rest / (real) zeros;
                }
            }
        }
         P /= P.Sum();
    }
#endif


#if 0
    // aka: I-BVMM -- best for many outcomes
    real S = alpha.Sum(); // = N_obs
    real N = 0; // N is the number of symbols
    for (int i=0; i<n_outcomes; ++i) {
        if (alpha(i)) {
            N += 1;
        }
    }
    real Z = (1 + N) * prior_alpha + S; // total dirichlet mass
    P = (alpha + prior_alpha) / Z;
    real n_zero_outcomes = n_outcomes - N;
    if (n_zero_outcomes > 0) {
        real SA = 1.0 / n_zero_outcomes;
        for (int i=0; i<n_outcomes; ++i) {
            if (alpha(i)==0) {
                P(i) *= SA;
            }
        }
    }
#endif
	alpha[y]++;
    // P(y | B_k) = P(y | B_k, h_k) P(h_k | B_k) + (1 - P(h_k | B_k)) P(y | B_{k-1})
    w = exp(log_w_prior + log_w); 


    total_probability = P[y] * w + (1 - w) * probability;
#if 0
    std::cout << depth << ": P(y|h_k)=" << P[y] 
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif
    //real posterior = w * P[y] / total_probability; // real posterior
    //real posterior = w; // fake posterior
    //real log_posterior = log(w) + log(P[y]) - log(total_probability);
    //log_w = log(posterior) - log_w_prior;
    log_w = log(w * P[y] /total_probability) - log_w_prior;

    // This sometimes doesn't work
    //log_w = log_w + log(P[y]) - log(total_probability);

    // Make sure we have enough observations to justify adding a
    // node. This means at least as many as total outcomes.
    //real threshold = (real) n_outcomes; 

    // Go deeper when there has been at least one observations
    // node. 
    //real threshold = sqrt((real) depth);//2;
    real threshold = log((real) depth);//2;

    // Always go deepr, no matter what
    //real threshold = 0; 

    // Go deeper if the context is long enough and the number of
    // observations justifies it.
	if (x != history.end() && S >  threshold) {
		int k = *x;
		++x;
		if (!next[k]) {
			next[k] = new Node(this);
		}
		total_probability = next[k]->Observe(history, x, y, total_probability);
	}

    N_obs++;

	return total_probability;
}


void ContextTree::Node::Show()
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

int ContextTree::Node::NChildren()
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

ContextTree::ContextTree(int n_branches_, int n_symbols_, int max_depth_)
	: n_branches(n_branches_),
	  n_symbols(n_symbols_),
	  max_depth(max_depth_),
	  history(max_depth)
{
	root = new Node(n_branches, n_symbols);
}

ContextTree::~ContextTree()
{
	delete root;
}

real ContextTree::Observe(int x, int y)
{
    history.push_back(x);
	return root->Observe(history, history.begin(), y, 0);
}

void ContextTree::Show()
{
    root->Show();
	std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTree::NChildren()
{
	return root->NChildren();
}
