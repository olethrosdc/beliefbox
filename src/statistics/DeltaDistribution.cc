/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DeltaDistribution.h"

DeltaUniformDistribution::DeltaUniformDistribution(const Vector& lower_bound_,
                                                   const Vector& upper_bound_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      mean((lower_bound + upper_bound)*0.5),
      T(0)
{
    n_dim = lower_bound.Size();
    assert(n_dim == upper_bound.Size());

}

real DeltaUniformDistribution::pdf(const Vector& x) const
{
    assert(x.Size() == n_dim);

    if (T == 0) {
        return 1.0/Volume(upper_bound - lower_bound);
    }
    if (x < lower_bound || x > upper_bound) {
        return 0;
    } else {
        if (L1Norm(&x,  &mean) < 1e-6) {
            return 1;
        } else {
            return 0;
        }
    }
}

real DeltaUniformDistribution::Observe(const Vector& x) 
{
    real p = pdf(x);
    if (T == 0) {
        assert(p > 0);
        mean = x;
    }
    T++;
    return p;
}
