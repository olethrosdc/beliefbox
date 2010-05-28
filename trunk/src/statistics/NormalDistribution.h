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
