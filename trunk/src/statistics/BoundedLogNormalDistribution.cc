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

#include "BoundedLogNormal.h"

BoundedLogNormal::BoundedLogNormal(const Vector& lower_bound, const Vector& upper_bound)
    : n_dim(lower_bound.Size()),
      a(lower_bound),
      b(upper_bound),
      c(0.5*(a+b)),
      d(0.5*(b-a))
{
    assert(a.Size() == b.Size());
    for (int i=0; i<n_dim; ++i) {
        assert(a(i) < b(i));
    }
    normal_density = new MultivariateNormalUnknownMeanPrecision(Vector(n_dim), 1.0, 1.0, Matrix::Unity(y_dim, y_dim));
        
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
/// Note that this the marginal likelihood!
real BoundedLogNormal::pdf(const Vector& x) const
{
    return normal_density->pdf(trasform(x));
}
/// The marginal log-likelihood
real BoundedLogNormal::logPdf(const Vector& x) const
{
    return normal_density->logPdf(trasform(x));
}
const Vector& BoundedLogNormal::getMean() const
{
    return inverse_transform(normal_density->getMean());
}

void BoundedLogNormal::calculatePosterior(const Vector& x)
{
    normal_density->calculatePosterior(x);
}
