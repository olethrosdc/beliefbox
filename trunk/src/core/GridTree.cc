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

#include "GridTree.h"

/** Construct a new grid.

    The grid is merely a mid-point subdivision of the n-dimensional space.
 */
GridTree::GridTree(Vector& lower_bound, Vector& upper_bound)
{
    assert(lower_bound.Size() == upper_bound.Size());
    root = new Node(lower_bound, upper_bound, 0);
}

/** Get the index of the interval containing x.
    
    For any \f$l, u \in R^n\f$, we define the function \f$f(\cdot | l,
    u) : R^n \to \{0, \ldots, 2^{n} -1\}\f$:
    \f[
    f(x | l, u) = \sum_{k=0}^{d - 1} 2^k I\left\{x_k > \frac{l_k + u_k}{2}\right\},
    \f]
    where \f$I\f$ is an indicator function.
 */
std::vector<int> GridTree::getInterval(Vector& x)
{
    std::vector<int> tmp;
    return tmp;
}
