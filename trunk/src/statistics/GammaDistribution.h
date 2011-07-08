/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of Sthe License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GAMMA_DISTRIBUTION_H
#define GAMMA_DISTRIBUTION_H

#include "Distribution.h"
#include "ExponentialDistribution.h"

/** Gamma distribution.
 *
 * The gamma distribution is conjugate with the exponential, Poisson
 * and the normal distribution for known mean.
 *
 * \f$\alpha, \beta\f$ are the shape and inverse scale parameters respectively.
 */
class GammaDistribution : public Distribution
{
public:
    real alpha;
    real beta;
    /// Beta distribution default constructor: alpha,beta=1, the
    /// uniform distribution
    GammaDistribution();
    GammaDistribution(real alpha_, real beta_);
    virtual real pdf(real x) const;
    virtual real log_pdf(real x) const;
    virtual real generate();
    virtual real generate() const;
    real setMaximumLikelihoodParameters(std::vector<real>& x, int n_iterations);
    real log_pdf(std::vector<real>& x) const;
};

/// Conjugate prior to gamma distribution with unknown shape and
/// inverse scale
class GammaDistributionUnknownShapeScale : public ConjugatePrior
{
public:
    real lambda; ///< for \f$P(\alpha | \lambda)$
    real mu; ///< for $P(\beta \mid \alpha, \mu, \nu)$
    real nu; ///< $P(\beta \mid \alpha, \mu, \nu)$
    real S; ///< sum of observations
    int T; ///< number of observations
    GammaDistributionUnknownShapeScale(real lambda_ = 1, real mu_ = 1, int nu_ = 1); ///< constructor
    virtual ~GammaDistributionUnknownShapeScale(); ///< destructor
    virtual void calculatePosterior(real x); ///< calculate posterior from an observation
    virtual real Observe(real x); ///< observe, return posterior..
    virtual real pdf (real x) const; ///< pdf
    /// observe sequence, return likelihood
    /// A number of samples must be specified
    virtual real Likelihood(std::vector<real>& x, int K) const
    {
        return exp(LogLikelihood(x, K));
    }
    /// observe sequence, return log likelihood
    virtual real LogLikelihood(std::vector<real>& x, int K) const; 
    virtual real generate() const;
};

#endif
