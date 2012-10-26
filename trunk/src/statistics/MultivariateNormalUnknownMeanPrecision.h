/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef MULTIVARIATE_NORMAL_UNKNOWN_MEAN_PRECISION_H
#define MULTIVARIATE_NORMAL_UNKNOWN_MEAN_PRECISION_H

#include "Student.h"

/** A multivariate normal distribution with unknown mean and precision.
    
    We have a sample \f$x_1, \ldots, x_n\f$ from a normal distribution
    with unknown mean \f$m\f$ and precision matrix \f$r\f$. We denote the
    normal density by \f$f(x \mid m, r)\f$. The prior joint parameter
    distribution \f$\xi_0(m, r) = \xi_0(m \mid r) \xi_0(r)\f$ is specified
    as follows.  The conditional distribution of \f$m\f$ given \f$r\f$
    is \f$\xi_0(m \mid r) = f(m \mid \mu_0, \tau_0 r)\f$ and the marginal of
    the precision is \f$\xi_0(r) = g(r \mid \alpha_0, T_0)\f$. 

    The predictive posterior distribution \f$\xi_n(x_{n+1})\f$, as
    well as the marginal mean distributions are actually both
    student-t distributions.
 */
class MultivariateNormalUnknownMeanPrecision : public VectorDistribution
{
protected:
    int n_dim; ///< number of dimensions
    Student marginal_mean; ///< the marginal distribution of m
    Student marginal; ///< the marginal distribution of x
public:
    // Parameters for m | R = r ~ Normal(\mu, \tau r)
    // where r is a given precision matrix
    Vector mu_0; ///< prior mean
    real tau_0; ///< prior accuracy
    Vector mu_n; ///< current mean
    real tau_n; ///< current accuracy

    // parameters for \xi(r) = Wishart(r | \alpha, T)
    real alpha_0; ///< prior alpha
    Matrix T_0; ///< prior T
    real alpha_n; ///< posterior alpha
    Matrix T_n; ///< posterior T

    // auxilliary parameters
    Vector bx_n; ///< \f$\bar{x}_n = \frac{1}{n} \sum_{i=1}^n x_i\f$.
    Matrix M_2n; ///< \f$M_{2,n} = \sum_{i=1}^n (x_n - \bar{x})^2\f$.
    int n;
    real sum;
    MultivariateNormalUnknownMeanPrecision();
    MultivariateNormalUnknownMeanPrecision(const Vector& mu, const real tau, const real alpha, const Matrix& T);
    void Reset();
    virtual ~MultivariateNormalUnknownMeanPrecision();
    virtual void generate(Vector& x) const
    {
        x = generate();
`    }
    virtual Vector generate();
    virtual Vector generate() const;
    /// This is the probability of a particular mean (ignoring the covariance)
    virtual real pdf(const Vector& mean) const; 
    /// Note that this the marginal likelihood!
    virtual real marginal_pdf(const Vector& x) const;
    /// The marginal log-likelihood
    virtual real log_pdf(const Vector& mean) const;
    virtual const Vector& getMean() const;
    virtual void calculatePosterior(const Vector& x);
    virtual real Observe(const Vector& x);
    void Show();
};

#endif
