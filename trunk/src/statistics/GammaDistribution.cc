/* -*- Mode: C++; -*- */
// VER: $Id: Distribution.c,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ranlib.h"
#include "GammaDistribution.h"
#include "SpecialFunctions.h"
#include "Distribution.h"
#include "ExponentialDistribution.h"

GammaDistribution::GammaDistribution() : alpha(1.0), beta(1.0)
{
}


GammaDistribution::GammaDistribution(real alpha_, real beta_) : alpha(alpha_), beta(beta_)
{
}

/** The log-pdf

    \f[
    \ln[ x^{\alpha -1} e^{-\beta x} \beta^\alpha / \Gamma(\alpha)]
    =
    (\alpha - 1) \ln x - \beta x + \alpha \ln \beta - \ln Gamma(\alpha)
    \f]
 */
real GammaDistribution::log_pdf(real x) const
{
   if (x<0.0 || alpha<0) {
       return LOG_ZERO;
    }
    real log_pdf = log(x)*(alpha-1.0) + log(beta)*alpha - beta*x - logGamma(alpha);
    return log_pdf;
}

/// Gamma has a support \f$[0, \infty)\f$, with pdf
/// \f$x^{\alpha-1} \frac{\beta^\alpha \exp(-\beta x)}{\Gamma(\alpha)}\f$,
/// with \f$k \geq 0\f$
real GammaDistribution::pdf(real x) const
{
    if (x<0.0 || alpha<0) {
        return 0.0;
    }
    real log_pdf = log(x)*(alpha-1.0) - log(beta)*alpha - beta*x - logGamma(alpha);
    return exp(log_pdf);
}

/// Generate using ranlib
real GammaDistribution::generate()
{
    return gengam(alpha, beta);
}

/// Generate using ranlib
real GammaDistribution::generate() const
{
    return gengam(alpha, beta);
}

/// Set the maximum likelihood parameters. Return the likelihood at that point.
///
/// Unfortunately this can only be done approximately. We generate
/// \f$\alpha\f$ randomly and subsequently find the ML value for
/// \f$\beta\f$ given \f$\alpha\f$. After a number of iterations, we
/// report the highest values.
real GammaDistribution::setMaximumLikelihoodParameters(const std::vector<real>& x, int n_iterations)
{
    real max_alpha = alpha;
    real max_beta = beta;
    real max_likelihood = Distribution::log_pdf(x);
    ExponentialDistribution Exp;
    real Z = ((real) x.size() / Sum(x));
    for (int iter=0; iter<n_iterations; ++iter) {
        alpha = Exp.generate();
        beta = alpha * Z;
        real likelihood = Distribution::log_pdf(x);
        if (likelihood > max_likelihood) {
            max_likelihood = likelihood;
            max_alpha = alpha;
            max_beta = beta;
        }
    }
    alpha = max_alpha;
    beta = max_beta;
    return max_likelihood;
}



// ----- conjugate  prior for gamma distribution ------ //

GammaDistributionUnknownShapeScale::GammaDistributionUnknownShapeScale(real lambda_, real mu_, int nu_) :
    lambda(lambda_),
    mu(mu_),
    nu(nu_),
    S(0),
    T(0)
{
}

GammaDistributionUnknownShapeScale::~GammaDistributionUnknownShapeScale()
{
}

void GammaDistributionUnknownShapeScale::calculatePosterior(real x)
{
    T ++;
    S += x;
}
real GammaDistributionUnknownShapeScale::Observe(real x) 
{
    real p = pdf(x);
    calculatePosterior(x);
    return p;
}
real GammaDistributionUnknownShapeScale::pdf (real x) const
{
    return 0.0;
}

real GammaDistributionUnknownShapeScale::LogLikelihood(std::vector<real>& x, int K) const
{
    real t = (real) T;
    ExponentialDistribution Exp(lambda);
    real log_likelihood = LOG_ZERO;
    for (int k=0; k<K; ++k) {
        real alpha = Exp.generate();
        //GammaDistribution prior_gamma(mu + t * alpha, nu + S);
        GammaDistribution prior_gamma(mu, nu);
        real beta = prior_gamma.generate();
        real log_p = 0.0;
        int n = x.size();
        GammaDistribution marginal_gamma(alpha, beta);
        for (int i=0; i<n; ++i) {
            log_p += marginal_gamma.log_pdf(x[i]);
        }
        log_likelihood = logAdd(log_likelihood, log_p);
    }
    log_likelihood -= log(K);
    return log_likelihood;
}

real GammaDistributionUnknownShapeScale::generate() const
{
    ExponentialDistribution Exp(lambda);
    real alpha = Exp.generate();
    GammaDistribution prior_gamma(mu + (real) T * alpha, nu + S);
    real beta = prior_gamma.generate();
    real log_p = 0.0;
    GammaDistribution marginal_gamma(alpha, beta);
    return marginal_gamma.generate();
}

