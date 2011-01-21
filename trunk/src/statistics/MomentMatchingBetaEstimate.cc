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

#include "MomentMatchingBetaEstimate.h"

MomentMatchingBetaEstimate::MomentMatchingBetaEstimate(const Vector& lower_bound, const Vector& upper_bound)
    : n_dim(lower_bound.Size()),
      T(0),
      a(lower_bound),
      b(upper_bound),
      d(b - a),
      c(Vector::Unity(n_dim) / d)
{
    assert(a.Size() == b.Size());
    for (int i=0; i<n_dim; ++i) {
        assert(a(i) < b(i));
        beta.push_back(new BetaDistribution(0.5, 0.5));
    }
}


void MomentMatchingBetaEstimate::Reset()
{
    for (int i=0; i<n_dim; ++i) {
        beta[i]->alpha = 0.5;
        beta[i]->beta = 0.5;
    }

}

MomentMatchingBetaEstimate::~MomentMatchingBetaEstimate()
{
    for (int i=0; i<n_dim; ++i) {
        delete beta[i];
    }
}

Vector MomentMatchingBetaEstimate::generate()
{
    Serror("Fix me!\n");
    return c;
} 

Vector MomentMatchingBetaEstimate::generate() const
{
    Serror("Fix me!\n");
    return c;
}

real MomentMatchingBetaEstimate::Observe(const Vector& x)
{
    real p = 1;
    Vector y = transform(x);
    for (int i=0; i<n_dim; ++i) {
        p *= beta[i]->pdf(y(i));
        
        real& a = beta[i]->alpha;
        real& b = beta[i]->beta;
        
    }
    return p;
}


/// Note that this the marginal likelihood!
real MomentMatchingBetaEstimate::pdf(const Vector& x) const
{
    real p = 1;
    Vector y = transform(x);
    for (int i=0; i<n_dim; ++i) {
        p *= beta[i]->pdf(y(i));
    }
    return p;
}

/// The marginal log-likelihood
real MomentMatchingBetaEstimate::log_pdf(const Vector& x) const
{
    real log_p = 0;
    Vector y = transform(x);
    for (int i=0; i<n_dim; ++i) {
        log_p += beta[i]->log_pdf(y(i));
    }
    return log_p;

}

const Vector& MomentMatchingBetaEstimate::getMean() const
{
    Serror("Fix me!\n");
    return Vector(n_dim);
}


Vector MomentMatchingBetaEstimate::transform(const Vector& x) const
{
    return (x - a) * c;
}

Vector MomentMatchingBetaEstimate::inverse_transform(const Vector& x) const
{
    return a + x * d;
}
