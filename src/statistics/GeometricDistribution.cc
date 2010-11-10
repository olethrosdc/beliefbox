/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "GeometricDistribution.h"
#include "Random.h"

real GeometricDistribution::generate()
{
    real x = urandom();
    return floor(log(x) * gen_scale);
}

real GeometricDistribution::pdf(real x) const
{
    if (x == floor(x)) {
        return p * pow((1 - p), x);
    }
    return 0.0;
}

