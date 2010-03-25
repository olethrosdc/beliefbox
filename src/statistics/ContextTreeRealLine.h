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

#ifndef CONTEXT_TREE_REAL_LINE_H
#define CONTEXT_TREE_REAL_LINE_H

#include <vector>
#include "real.h"
#include "Vector.h"
#include "Ring.h"


/** Context tree on the real line.

    This is a recursive version of the structure described in Marcus
    Hutter's Bayesian Infinite Tree paper.

    The BVMM approach (implemented in ContextTree.[h|cc]) is not applicable
    here.
*/
class ContextTreeRealLine
{
public:
	// public classes
	struct Node
	{
        real lower_bound; ///< looks at x > lower_bound
        real upper_bound; ///< looks at x < upper_bound
        real new_bound; ///< how to split
        const int n_branches; ///< number of branches
        const int depth; ///< depth of the node
        const int max_depth; ///< maximum depth
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
		real P; ///< probability of next symbols
        Vector alpha; ///< number of times seen in each quadrant
        real w; ///< backoff weight
        real log_w; ///< log of w
        real log_w_prior; ///< initial value
		Node(real lower_bound_,
             real upper_bound_,
             int n_branches_,
             int max_depth_);
		Node(Node* prev_, real lower_bound_, real upper_bound_);
		~Node();
		real Observe(real x,
					 real probability);
		real pdf(real x,
                 real probability);
		void Show();
		int NChildren();	
        int S;
	};
	
	// public methods
	ContextTreeRealLine(int n_branches_ = 2, int max_depth_= 0, real lower_bound = 0, real upper_bound = 1);
	~ContextTreeRealLine();
	real Observe(real x);
	real pdf(real x);
	void Show();
	int NChildren();
protected: 
	int n_branches;
	int max_depth;
	Node* root;
};



#endif
