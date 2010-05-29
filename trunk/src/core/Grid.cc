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

#include "Grid.h"

/** Construct a new grid.

    The grid is merely a mid-point subdivision of the n-dimensional space.
 */
Grid::Grid(Vector& lower_bound_, Vector& upper_bound_)
    :
    lower_bound(lower_bound_),
    upper_bound(upper_bound_)
{
    assert(lower_bound.Size() == upper_bound.Size());
    mid_point = (lower_bound + upper_bound)*0.5;
    n_dimensions = lower_bound.Size();
}

/** Get the index of the interval containing x.
    
    For any \f$l, u \in R^n\f$, we define the function \f$f(\cdot | l,
    u) : R^n \to \{0, \ldots, 2^{n} -1\}\f$:
    \f[
    f(x | l, u) = \sum_{k=0}^{d - 1} 2^k I\left\{x_k > \frac{l_k + u_k}{2}\right\},
    \f]
    where \f$I\f$ is an indicator function.
 */
int Grid::getInterval(Vector& x)
{
    int d = 1;
    int y = 0;
    for (int i=0; i<n_dimensions; ++i) {
        assert(x[i] >= lower_bound[i]);
        assert(x[i] <= upper_bound[i]);        
        if (x[i] > mid_point[i]) {
            y += d;
        }
        d <<= 1;
    }
    return y;
}



/** Construct a new grid.

    Construct a new grid, subdividing each dimension in K parts.
 */
EvenGrid::EvenGrid(Vector& lower_bound_, Vector& upper_bound_, int K_)
    :
    lower_bound(lower_bound_),
    upper_bound(upper_bound_),
    K(K_)
{
    assert(lower_bound.Size() == upper_bound.Size());
    assert(K >= 1);
    delta = (upper_bound - lower_bound) / K ;
    n_dimensions = lower_bound.Size();
    n_intervals = (int) floor(pow(K, n_dimensions));
}

/** Get the index of the interval containing x.
    
    For any \f$l, u \in R^n\f$ and \f$k \in \{1, 2, \ldots \}\f$, we define the function \f$f(\cdot | k, l, u) : R^n \to \{0, \ldots, k^{n} -1\}\f$:
    \f[
    f(x | k, l, u) = \sum_{d=0}^n k^d I_d(x | k, l, u)
    \f]
    where \f$I(\cdot | k, l, u) R^n \to \{0, \ldots, k - 1\}^d \f$ 
    is an indicator function such that \f$I_d(x | k, l, u) = \max\{0, \min[k - 1, s_d]\}\f$ and \f$s_d = \lfloor \frac{k (x_d - l_d)}{u_d - l_d} \rfloor \f$.


    Note that each one of the extreme intervals implicitly includes
    plus or minus infinity. 
 */
int EvenGrid::getInterval(Vector& x) const
{
    int d = 1;
    int y = 0;

    for (int i=0; i<n_dimensions; ++i) {
        real dx = (x[i] - lower_bound[i]) / delta[i];
        int n = (int) floor(dx);
        if (n >= K) {
            n = K - 1;
        } else if (n < 0) {
            n = 0;
        }
        y += d * n;
        d *= K;
    }
    assert(y >= 0 && y < n_intervals);
    return y;
}

/// Get the number of intervals.
int EvenGrid::getNIntervals() const
{
    return n_intervals;
}
