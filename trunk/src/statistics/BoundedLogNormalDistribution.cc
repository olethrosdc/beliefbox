/* -*- Mode: C++ -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BoundedLogNormalDistribution.h"

BoundedLogNormal::BoundedLogNormal(const Vector& lower_bound, const Vector& upper_bound)
    : n_dim(lower_bound.Size()),
      a(lower_bound),
      b(upper_bound),
      c((a+b) * 0.5),
      d(pow(b-a, -2.0))
{
    assert(a.Size() == b.Size());
    for (int i=0; i<n_dim; ++i) {
        assert(a(i) < b(i));
    }
    normal_density =
        new MultivariateNormalUnknownMeanPrecision(Vector(n_dim), 
                                                   1.0, 
                                                   1.0,
                                                   Matrix::Unity(n_dim, n_dim));
        
}


void BoundedLogNormal::Reset()
{
    normal_density->Reset();
}

BoundedLogNormal::~BoundedLogNormal()
{
    delete normal_density;
}
Vector BoundedLogNormal::generate()
{
    Serror("Fix me!\n");
    return c;
} 
Vector BoundedLogNormal::generate() const
{
    Serror("Fix me!\n");
    return c;
}

real BoundedLogNormal::Observe(const Vector& x)
{
    return normal_density->Observe(transform(x));
}


/// Note that this the marginal likelihood!
real BoundedLogNormal::pdf(const Vector& x) const
{
    return normal_density->pdf(transform(x));
}
/// The marginal log-likelihood
real BoundedLogNormal::logPdf(const Vector& x) const
{
    return normal_density->logPdf(transform(x));
}
const Vector& BoundedLogNormal::getMean() const
{
    Serror("Fix me!\n");
    return normal_density->getMean();
}

void BoundedLogNormal::calculatePosterior(const Vector& x)
{
    normal_density->calculatePosterior(x);
}
Vector BoundedLogNormal::transform(const Vector& x) const
{
    return (c - x) * log (d * (x - a) * (b - x));
}
Vector BoundedLogNormal::inverse_transform(const Vector& x) const
{
    Serror("Fix me!\n");
    return x;
}
