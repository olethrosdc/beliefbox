/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "ExponentialDistribution.h"
#include "Distribution.h"


real ExponentialDistribution::generate() const
{
    real x = urandom();
    return - log (1.0 - x) / l;
}

real ExponentialDistribution::pdf(real x) const
{
    real d = x - m;
    if (d>0.0) {
        return l * exp (-l*d);
    }
    return 0.0;
}

real ExponentialDistribution::log_pdf(real x) const
{
    real d = x - m;
    if (d>0.0) {
        return log(l) - l*d;
    }
    return LOG_ZERO;
}
