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
#ifndef MULTIVARIATE_NORMAL_H
#define MULTIVARIATE_NORMAL_H

#include "NormalDistribution.h"

/// Multivariate Gaussian probability distribution
class MultivariateNormal : public VectorDistribution
{
 private:
    int n_dim;
    Vector mean;
    Matrix accuracy;
    real determinant;
 public:
    MultivariateNormal(const int n_dim_);
    MultivariateNormal(const Vector& mean_, const Matrix& accuracy_);
    void setMean(const Vector& mean_)
    {
        mean = mean_;
    }
    void setAccuracy(const Matrix& accuracy_)
    {
        accuracy = accuracy_;
        accuracy.LUDecomposition(determinant);
    }
    virtual ~MultivariateNormal() {}
    virtual void generate(Vector& x) const;
    virtual Vector generate() const;
    virtual real log_pdf(const Vector& x) const;
    virtual real pdf(const Vector& x) const
    {
        return exp(log_pdf(x));
    }
    void Show() const;
};


/** A multivariate normal distribution with unknown mean and precision.
    
    We have a sample \f$x_1, \ldots, x_n\f$ from a normal distribution
    with unknown mean \f$m\f$ and precision matrix \f$r\f$. We denote the
    normal density by \f$f(x \mid m, r)\f$. The prior joint parameter
    distribution \f$\xi_0(m, r) = \xi_0(m \mid r) \xi_0(r)\f$ is specified
    as follows.  The conditional distribution of \f$m\f$ given \f$r\f$
    is \f$\xi_0(m \mid r) = f(m \mid \mu_0, \tau_0 r)\f$ and the marginal of
    the precision is \f$\xi_0(r) = g(r \mid \alpha_0, T_0)\f$. 

    The predictive posterior distribution \f$\xi_n(x_{n+1})\f$ is
    actually a generalised student-t distribution, but here we are
    hacking it as a normal.
 */
class MultivariateNormalUnknownMeanPrecision
{
protected:
    int n_dim;
    MultivariateNormal p_x_mr; ///< \f$\xi(x) = \int f(x | m, r) \, d\xi(m ,r)\f$  
public:
    // paramters for \xi(m | r) = f(m | \mu, \tau r)
    Vector mu_0; ///< prior mean
    real tau_0; ///< prior accuracy
    Vector mu_n; ///< current mean
    real tau_n; ///< current accuracy

    // parameters for \xi(r) = g(r | \alpha, \beta)
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
    virtual Vector generate();
    virtual Vector generate() const;
    /// Note that this the marginal likelihood!
    virtual real pdf(const Vector& x) const;
    /// The marginal log-likelihood
    virtual real log_pdf(const Vector& x) const;
    virtual const Vector& getMean() const;
    virtual void calculatePosterior(const Vector& x);
    virtual real Observe(const Vector& x);
    void Show();
};

#endif
