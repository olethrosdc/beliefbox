/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef NORMAL_DISTRIBUTION_H
#define NORMAL_DISTRIBUTION_H

#include "Distribution.h"
#include "Matrix.h"
#include "Vector.h"


/// Gaussian probability distribution 
class NormalDistribution : public ParametricDistribution {
 private:
    mutable bool cache;
    mutable real normal_x, normal_y, normal_rho;
 public:
    real m; ///< mean
    real s; ///< standard deviation
    NormalDistribution() {m=0.0; s=1.0; cache = false;}
    /// Normal dist. with given mean and std
    NormalDistribution(real mean, real std)
    {
        setMean (mean);
        setVariance (std*std);
    }
    virtual Distribution* clone ()
    {
        NormalDistribution* d  = new NormalDistribution;
        d->m = m;
        d->s = s;
        return d;
    }
    virtual ~NormalDistribution() {}
    virtual real generate();
    virtual real generate() const;
    virtual real pdf(real x) const;
    void setSTD(real std)
    {
        s = std;
    }
    virtual void setVariance (real var) 
    {
        s = sqrt(var);
    } 
    virtual void setMean (real mean)
    {
        m = mean;
    }
    virtual real getMean () const
    {
        return m;
    }
};


class NormalDistributionUnknownMean : public ConjugatePrior
{
protected:
    NormalDistribution observations;
    NormalDistribution prior;
public:
    real mu_0; ///< prior mean
    real tau_0; ///< prior accuracy
    real mu_n; ///< current mean
    real tau_n; ///< current accuracy
    real tau; ///< observation accuracy
    int n;
    real sum;
    NormalDistributionUnknownMean() {
        mu_n = 0.0;
        tau_n =1.0;
        tau = 1.0;
        Reset();
    }

    NormalDistributionUnknownMean(real mu_0_, real tau_0_, real tau_)
        : mu_0(mu_0_), tau_0(tau_0_), tau(tau_)
    {
        Reset();
    }
    void Reset()
    {
        observations.setMean(mu_0);
        observations.setSTD(1.0 / tau);
        prior.setMean(mu_0);
        prior.setVariance(1.0 / (tau_0 * tau_0));
        n = 0;
        sum = 0.0;
        tau_n = tau_0;
        mu_n = mu_0 * tau_0;
    }
    virtual ~NormalDistributionUnknownMean()
    {
    }
    virtual real generate();
    virtual real generate() const;
    /// Note that this the marginal likelihood!
    virtual real pdf(real x) const;
    virtual real getMean() const;
    virtual void calculatePosterior(real x);
    real Observe(real x)
    {
        real p = pdf(x);
        calculatePosterior(x);
        return p;
    }
};

/** A normal distribution with unknown mean and precision.
    
    We have a sample \f$x_1, \ldots, x_n\f$ from a normal distribution
    with unknown mean \f$m\f$ and precision \f$r\f$. We denote the
    normal density by \f$f(x \mid m, r)\f$. The prior joint parameter
    distribution \f$\xi_0(m, r) = \xi_0(m \mid r) \xi_0(r)\f$ is specified
    as follows.  The conditional distribution of \f$m\f$ given \f$r\f$
    is \f$\xi_0(m \mid r) = f(m \mid \mu_0, \tau_0 r)\f$ and the marginal of
    the precision is \f$\xi_0(r) = g(r \mid \alpha_0, \beta_0\f$. 
 */
class NormalUnknownMeanPrecision: public ConjugatePrior
{
protected:
    NormalDistribution p_x_mr; ///< \f$f(x | m, r)\f$
    //NormalDistribution p_m_r; ///< \f$\xi(m | r) = f(m \mid \mu, \tau r)\f$
    //GammaDistribution p_r; ///< \f$\xi(r) = g(r \mid \alpha, \beta)\f$
public:
    // paramters for \xi(m | r) = f(m | \mu, \tau r)
    real mu_0; ///< prior mean
    real tau_0; ///< prior accuracy
    real mu_n; ///< current mean
    real tau_n; ///< current accuracy

    // parameters for \xi(r) = g(r | \alpha, \beta)
    real alpha_0; ///< prior alpha
    real beta_0; ///< prior beta
    real alpha_n; ///< posterior alpha
    real beta_n; ///< posterior beta

    // auxilliary parameters
    real bx_n; ///< \f$\bar{x}_n = \frac{1}{n} \sum_{i=1}^n x_i.
    real M_2n; ///< \f$M_{2,n} = \sum_{i=1}^n (x_n - \bar{x})^2\f$.
    int n;
    real sum;
    NormalUnknownMeanPrecision()
    {
        mu_n = 0.0;
        tau_n =1.0;
        Reset();
    }

    NormalUnknownMeanPrecision(real mu_0_, real tau_0_, real tau_)
        : mu_0(mu_0_), tau_0(tau_0_)
    {
        Reset();
    }
    void Reset()
    {
        p_x_mr.setMean(mu_0);
        p_x_mr.setVariance(1.0 / (tau_0 * tau_0));
        n = 0;
        sum = 0.0;
        tau_n = tau_0;
        mu_n = mu_0 * tau_0;
        M_2n = 0;
        bx_n = 0;
    }
    virtual ~NormalUnknownMeanPrecision()
    {
    }
    virtual real generate()
    {
        Serror("Fix me!\n");
        return mu_n;
    }
    virtual real generate() const
    {
        Serror("Fix me!\n");
    }
    /// Note that this the marginal likelihood!
    virtual real pdf(real x) const;
    virtual real getMean() const;
    virtual void calculatePosterior(real x);
    real Observe(real x)
    {
        real p = pdf(x);
        calculatePosterior(x);
        return p;
    }
};




/// Multivariate Gaussian probability distribution
class MultivariateNormal : public VectorDistribution {
 private:
    Vector mean;
    Matrix accuracy;
 public:

    MultivariateNormal(); 
    MultivariateNormal(Vector mean_, Matrix std_);
    virtual ~MultivariateNormal() {}
    virtual void generate(Vector& x);
    virtual Vector generate();
    virtual real pdf(Vector& x) const;
};

#endif
