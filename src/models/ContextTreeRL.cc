
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

#include "ContextTreeRL.h"

ContextTreeRL::Node::Node(int n_branches_,
                          int n_outcomes_)
    : n_branches(n_branches_),
      n_outcomes(n_outcomes_),
      depth(0),
      prev(NULL),
      next(n_branches),
      P(n_outcomes), alpha(n_outcomes), prior_alpha(0.5),
      w(1), log_w(0), log_w_prior(0), Q(0)
{
    for (int i=0; i<n_outcomes; ++i) {
        P(i) = 1.0 / (real) n_outcomes;
        alpha(i) = 0;
    }
}

/// Make a node for K symbols at nominal depth d
ContextTreeRL::Node::Node(ContextTreeRL::Node* prev_)
    : n_branches(prev_->n_branches),
      n_outcomes(prev_->n_outcomes),
      depth(prev_->depth + 1),
      prev(prev_),
      next(n_branches),
      P(n_outcomes),
      alpha(n_outcomes),
      prior_alpha(0.5),
      log_w(0),
      log_w_prior(prev_->log_w_prior - log(2)),
      Q(0)
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
ContextTreeRL::Node::~Node()
{
    for (int i=0; i<n_branches; ++i) {
        delete next[i];
    }
}

real ContextTreeRL::Node::Observe(Ring<int>& history,
                                  Ring<int>::iterator x,
                                  int y,
                                  real r,
                                  real probability)
{
    real total_probability = 0;
    // calculate probabilities

    // Standard
#if 0
    real S = alpha.Sum();
    real Z = 1.0 / (prior_alpha * (real) n_outcomes + S);
    P = (alpha + prior_alpha) * Z;
    P /= P.Sum();
    //P[y] = (alpha[y] + prior_alpha) * Z;
#endif

#if 1
    // aka: I-BVMM -- best for many outcomes
    real S = alpha.Sum();
    real N = 0;
    for (int i=0; i<n_outcomes; ++i) {
        if (alpha(i)) {
            N += 1;
        }
    }
    real Z = (1 + N) * prior_alpha + S;
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

    // Do it for probability too
    real p_reward = reward_prior.Observe(r);
        
    // P(y | B_k) = P(y | B_k, h_k) P(h_k | B_k) + (1 - P(h_k | B_k)) P(y | B_{k-1})
    w = exp(log_w_prior + log_w); 
    
    
    real p_observations = P[y] * p_reward;
    total_probability = p_observations * w + (1 - w) * probability;
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
    log_w = log_w + log_w_prior + log(p_observations) - log(total_probability) - log_w_prior;

    // Make sure we have enough observations to justify adding a
    // node. This means at least as many as total outcomes.
    //real threshold = (real) n_outcomes; 

    // Go deeper when there has been at least one observations
    // node. 
    real threshold = 2;

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
        total_probability = next[k]->Observe(history, x, y, r, total_probability);
    }

    return total_probability;
}


void ContextTreeRL::Node::Show()
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

int ContextTreeRL::Node::NChildren()
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

ContextTreeRL::ContextTreeRL(int n_branches_, int n_symbols_, int max_depth_)
    : n_branches(n_branches_),
      n_symbols(n_symbols_),
      max_depth(max_depth_),
      history(max_depth)
{
    root = new Node(n_branches, n_symbols);
}

ContextTreeRL::~ContextTreeRL()
{
    delete root;
}

real ContextTreeRL::Observe(int x, int y, real r)
{
    history.push_back(x);
    return root->Observe(history, history.begin(), y, r, 0);
}

void ContextTreeRL::Show()
{
    root->Show();
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ContextTreeRL::NChildren()
{
    return root->NChildren();
}
void ContextTreeRL::QLearning(real step_size, int action, int observation, real reward)
{
    
}
