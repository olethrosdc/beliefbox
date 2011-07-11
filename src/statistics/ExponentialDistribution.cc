/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "ExponentialDistribution.h"
#include "Distribution.h"
#include "SpecialFunctions.h"


real ExponentialDistribution::generate() const
{
    real x = urandom();
    return - log (1.0 - x) / l;
}

real ExponentialDistribution::pdf(real x) const
{
    real d = x - m;
    if (d>0.0) {
        return l * exp (-l*d);
    }
    return 0.0;
}

real ExponentialDistribution::log_pdf(const std::vector<real>& x) const
{
    int T = x.size();
    real l = 0;
    for (int i=0; i<T; ++i) {
        l += log_pdf(x[i]);
    }
    return l;
}

real ExponentialDistribution::log_pdf(real x) const
{
    real d = x - m;
    if (d>0.0) {
        return log(l) - l*d;
    }
    return LOG_ZERO;
}


/// There is a closed form solution for this.. let's do it!
real ExponentialDistribution::setMaximumLikelihoodParameters(const std::vector<real>& x)
{
    setMean(Sum(x)  / (real) x.size());
}

// ----- gamma distribution conjugate prior for the exponential parameter ----- //



GammaExponentialPrior::GammaExponentialPrior(real alpha, real beta)
: gamma_prior(alpha, beta)
{
    GammaDistribution test;
}


real GammaExponentialPrior::Observe(real x)
{
    calculatePosterior(x);
    return 0;
}

void GammaExponentialPrior::calculatePosterior(real x)
{
    gamma_prior.alpha += 1.0;
    gamma_prior.beta += x;
}

real GammaExponentialPrior::log_pdf(real x, int n_iterations) const
{
    real log_likelihood = 0;
    for (int iter = 0; iter < n_iterations; ++iter) {
        real lambda = gamma_prior.generate();
        ExponentialDistribution exponential_sample(lambda);
        log_likelihood = logAdd(log_likelihood, exponential_sample.log_pdf(x));
    }
    log_likelihood -= log(n_iterations);

    real alpha = gamma_prior.alpha;
    real beta = gamma_prior.beta;

    real l2 = logGamma(alpha + 1) - logGamma(alpha) + alpha * log(beta) - (alpha + 1) * log(beta + x);
    printf ("%f %f\n", log_likelihood, l2);
    return log_likelihood;
}

real GammaExponentialPrior::generate() const
{
    real lambda = gamma_prior.generate();
    ExponentialDistribution exponential_sample(lambda);
    return exponential_sample.generate();
}


real GammaExponentialPrior::LogLikelihood(const std::vector<real>& x,
                                          int n_iterations) const
{
    real log_likelihood = 0;
    for (int iter = 0; iter < n_iterations; ++iter) {
        real lambda = gamma_prior.generate();
        ExponentialDistribution exponential_sample(lambda);
        log_likelihood = logAdd(log_likelihood, exponential_sample.log_pdf(x));
    }
    log_likelihood -= log(n_iterations);

    real alpha = gamma_prior.alpha;
    real beta = gamma_prior.beta;

    real l2 = 0;
    for (int i=0; i<(int) x.size(); ++i) {
        l2 += logGamma(alpha + 1) - logGamma(alpha) + alpha * log(beta) - (alpha + 1) * log(beta + x[i]);
    }
    //printf ("%f %f\n", log_likelihood, l2);
    return log_likelihood;
}
