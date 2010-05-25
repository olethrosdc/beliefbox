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

#ifndef GRID_TREE_H
#define GRID_TREE_H


#include "Grid.h"
#include <cassert>

/** A simple grid structure.

    It subdivides an \f$n\f$-dimensional space in \f$2^n\f$ subspace
    of equal volumes.
    
 */
struct GridTree
{
    struct Node
    {
        Grid grid;
        std::vector<GridTree::Node*> next;
        Node(Vector& lower_bound, Vector& upper_bound)
            : grid(lower_bound, upper_bound)
              next(lower_bound.Size())
        {

        }
    };
    Node* root;
    
    int n_dimensions;
    GridTree(Vector& lower_bound_, Vector& upper_bound_);
    std::vector<int> getInterval(Vector& x);
};



#endif
