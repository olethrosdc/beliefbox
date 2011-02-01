/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef POINT_PAIR
#define POINT_PAIR

#include "Vector.h"

/** A pair of two points. 

    This class is useful for conditional density estimation and
    regression.
 */
class PointPair
{
public:
    Vector x; ///< The first point
    Vector y; ///< The secoind point
    /// Constructor
    PointPair(const Vector& x_, const Vector& y_) : x(x_), y(y_)
    {
    }
};

#endif
