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

/** An alternative context tree implementation */
class ContextTree
{
public:

	// public classes

	struct Node
	{
		int depth; ///< depth
		std::vector<Node*> next; ///< pointers to next nodes
		Vector P; ///< probability of next symbols
		real w;
		/// Make a node for K symbols at nominal depth d
		Node(int d, int K) : depth(d), next(K), P(K)
		{
			real w = 0;
			for (int i=0; i<K; ++i) {
				next[i] = NULL;
			}
		}
	};
	
	// public methods

	ContextTree(int n_symbols_, int max_depth_= 0);
	virtual ~ContextTree();
	real Observe(int x);

protected: 
	int n_symbols;
	int max_depth;
	Node* root;
};


#endif
