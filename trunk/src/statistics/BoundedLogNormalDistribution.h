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

#ifndef BOUNDED_LOGNORMAL_DISTRIBUTION_H
#define BOUNDED_LOGNORMAL_DISTRIBUTION_H

#include "NormalDistribution.h"
#include "Matrix.h"
#include "Vector.h"
#include "Student.h"


/** A multivariate lognormal distribution with unknown mean and precision.
    
    We define the transform:
    \f[
    g_i(y) = (c_i - y_i)
    \ln \left[
    (d_i)^{-2} (y_i - a_i)(b_i - x_i)
    \right]
    \f] 
    where \f$c = \frac{a+b}{2}\f$ and \f$d = \frac{b-a}{2}\f$.
    Then \f$x_t = g(y_t)\f$, are normally distributed with unknown
    mean \f$m\f$ and precision matrix \f$r\f$. 
 */
class BoundedLogNormal
{
protected:
    MultivariateNormalUnknownMeanPrecision* normal_density;
    int n_dim;
    Vector transform(const Vector& x) const;
    Vector inverse_transform(const Vector& x) const;
public:
    Vector a; ///< lower bound
    Vector b; ///< upper bound
    Vector c; ///< center
    Vector d; ///< diff
    // auxilliary parameters
        
    BoundedLogNormal(const Vector& lower_bound, const Vector& upper_bound);
    void Reset();
    virtual ~BoundedLogNormal();
    virtual Vector generate();
    virtual Vector generate() const;
    /// Note that this the marginal likelihood!
    virtual real pdf(const Vector& x) const;
    /// The marginal log-likelihood
    virtual real logPdf(const Vector& x) const;
    virtual const Vector& getMean() const;
    virtual void calculatePosterior(const Vector& x);
    real Observe(const Vector& x);
    void Show();
};



#endif
