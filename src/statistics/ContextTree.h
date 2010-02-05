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

#ifndef CONTEXT_TREE_H
#define CONTEXT_TREE_H

#include <vector>
#include "real.h"
#include "Vector.h"
#include "Ring.h"


/** An alternative context tree implementation.

    Somewhat more generic.
    
    
 */
class ContextTree
{
public:
	// public classes
	struct Node
	{
        int n_branches;
        int n_outcomes;
		int depth; ///< depth
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
		Vector P; ///< probability of next symbols
		Vector alpha; ///< parameters of next symbols
        const real prior_alpha;
        real w;
        real log_w;
        real log_w_prior;
		Node(int n_branches_,
             int n_outcomes_)
            : n_branches(n_branches_),
              n_outcomes(n_outcomes_),
              depth(0),
              prev(NULL),
              next(n_branches),
              P(n_outcomes), alpha(n_outcomes), prior_alpha(0.5),
              w(0), log_w(LOG_ZERO), log_w_prior(0)
        {
            for (int i=0; i<n_outcomes; ++i) {
                P(i) = 1.0 / (real) n_outcomes;
                alpha(i) = 0;
            }
        }

		/// Make a node for K symbols at nominal depth d
		Node(Node* prev_)
            : n_branches(prev_->n_branches),
              n_outcomes(prev_->n_outcomes),
              depth(prev_->depth + 1),
              prev(prev_),
              next(n_branches),
              P(n_outcomes),
              alpha(n_outcomes),
              prior_alpha(0.5),
              w(0),
              log_w(LOG_ZERO),
              log_w_prior(prev_->log_w_prior - 2)
		{
			for (int i=0; i<n_branches; ++i) {
				next[i] = NULL;
			}
            for (int i=0; i<n_outcomes; ++i) {
                P(i) = 1.0 / (real) n_outcomes;
                alpha(i) = 0;
            }

		}

        /// make sure to kill all
        ~Node()
        {
            for (int i=0; i<n_branches; ++i) {
                delete next[i];
            }
        }
        real Observe(Ring<int>& history,
                Ring<int>::iterator x,
                int y,
                real probability)
        {
            real total_probability = 0;
            // calculate probabilities
            real Z = 1.0 / alpha.Sum();
            P = alpha * Z;

            ++x;
            if (x != history.end()) {
                if (!next[*x]) {
                    next[*x] = new Node(this);
                }
                total_probability = next[*x]->Observe(history, x, y, P[y]);
            }

            return total_probability;
        }
	};
	
	// public methods
	ContextTree(int n_symbols_, int max_depth_= 0);
	virtual ~ContextTree();
	real Observe(int x, int y);

protected: 
	int n_symbols;
	int max_depth;
	Node* root;
    Ring<int> history;
};


#endif
