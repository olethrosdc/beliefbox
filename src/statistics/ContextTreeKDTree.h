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

#ifndef CONTEXT_TREE_KD_TREE_H
#define CONTEXT_TREE_KD_TREE_H

#include <vector>
#include "real.h"
#include "Vector.h"
#include "Ring.h"


/** Context tree on the real line.

    This is a generalisation of the binary tree on [a,b] implemented in
	ContextTreeRealLine to a KD-tree.
*/
class ContextTreeKDTree
{
public:
    // public classes
    struct Node
    {
        Vector lower_bound; ///< looks at x > lower_bound
        Vector upper_bound; ///< looks at x < upper_bound
        real mid_point; ///< how to split
		int splitting_dimension; ///< dimension on which to do the split.

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

        Node(Vector& lower_bound_,
             Vector& upper_bound_,
             int n_branches_,
             int max_depth_);
        Node(Node* prev_, Vector& lower_bound_, Vector& upper_bound_);
        ~Node();
        real Observe(Vector& x,
                     real probability);
        real pdf(Vector& x,
                 real probability);
        void Show();
        int NChildren();    
        int S;
    };
    
    // public methods
    ContextTreeKDTree(int n_branches_, int max_depth_, Vector& lower_bound, Vector& upper_bound);
    ~ContextTreeKDTree();
    real Observe(Vector& x);
    real pdf(Vector& x);
    void Show();
    int NChildren();
protected: 
    int n_branches;
    int max_depth;
    Node* root;
};



#endif
