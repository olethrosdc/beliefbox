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

#ifndef GRID_H
#define GRID_H


#include "Vector.h"
#include <cassert>

/** A simple grid structure.

    It subdivides an \f$n\f$-dimensional space in \f$2^n\f$ subspaces.
 */
struct Grid
{
    Vector lower_bound;
    Vector upper_bound;
    Vector mid_point;
    int n_dimensions;
    Grid(Vector& lower_bound_, Vector& upper_bound_);
    int getInterval(Vector& x);
};

/** A simple grid structure.

    It subdivides an \f$n\f$-dimensional space in \f$K^n\f$ subspace.    
 */
struct EvenGrid
{
    Vector lower_bound;
    Vector upper_bound;
    Vector delta;
    int n_dimensions;
    int n_intervals;
    int K; ///< Number of subdivisions per dimension.
    EvenGrid(Vector& lower_bound_, Vector& upper_bound_, int K_);
    int getInterval(Vector& x) const;
    int getNIntervals() const;
 
};
#endif
