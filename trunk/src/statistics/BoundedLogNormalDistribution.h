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
    4(d_i)^{-2} (y_i - a_i)(b_i - x_i)
    \right]
    \f] 
    where \f$c = \frac{a+b}{2}\f$ and \f$d = \frac{b-a}{2}\f$.
    Then \f$x_t = g(y_t)\f$, are normally distributed with unknown
    mean \f$m\f$ and precision matrix \f$r\f$. We denote the normal
    density by \f$f(x \mid m, r)\f$. The prior joint parameter
    distribution \f$\xi_0(m, r) = \xi_0(m \mid r) \xi_0(r)\f$ is
    specified as follows.  The conditional distribution of \f$m\f$
    given \f$r\f$ is \f$\xi_0(m \mid r) = f(m \mid \mu_0, \tau_0 r)\f$
    and the marginal of the precision is \f$\xi_0(r) = g(r \mid
    \alpha_0, T_0)\f$.

    The predictive posterior distribution \f$\xi_n(x_{n+1})\f$ is
    actually a generalised student-t distribution, but here we are
    hacking it as a normal.
 */
class BoundedLogNormal
{
protected:
    MultivariateNormalUnknownMeanPrecision* normal_density;
    int n_dim;
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
