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

#include "NormalDistribution.h"
#include "ExponentialDistribution.h"
#include "Random.h"
#include "SpecialFunctions.h"

// Taken from numerical recipes in C
real NormalDistribution::generate() const
{
    if(!cache) {
        normal_x = urandom();
        normal_y = urandom();
        normal_rho = sqrt(-2.0 * log(1.0 - normal_y));
        cache = true;
    } else {
        cache = false;
    }
    
    if (cache) {
        return normal_rho * cos(2.0 * M_PI * normal_x) * s + m;
    } else {
        return normal_rho * sin(2.0 * M_PI * normal_x) * s + m; 
    }
}

real NormalDistribution::generate() 
{
    if(!cache) {
        normal_x = urandom();
        normal_y = urandom();
        normal_rho = sqrt(-2.0 * log(1.0 - normal_y));
        cache = true;
    } else {
        cache = false;
    }
    
    if (cache) {
        return normal_rho * cos(2.0 * M_PI * normal_x) * s + m;
    } else {
        return normal_rho * sin(2.0 * M_PI * normal_x) * s + m; 
    }
}
/// Normal distribution log-pdf
real NormalDistribution::log_pdf(real x) const
{
    real d = (m-x)/s;
    return -0.5 * d*d - 0.5 * log(2.0 * M_PI * s * s);
}

/// Normal distribution pdf
real NormalDistribution::pdf(real x) const
{
    real d = (m-x)/s;
    return exp(-0.5 * d*d)/(sqrt(2.0 * M_PI) * s);
}

/// Normal distribution Show
void NormalDistribution::Show() const
{
    printf("Normal - ");
    printf("Mean: %f Variance: %f, Accuracy: %f\n", m, s*s , 1.0/(s*s));
}


/// Set Maximum Likelihood parameters
///
/// Returns the log-likelihood of the parameters
real NormalDistribution::setMaximumLikelihoodParameters (std::vector<real>& x)
{
    int T = x.size();
    real t = (real) T;
    m = Sum(x) / t;
    s = 0;
    for (int i=0; i<T; ++i) {
        real d = x[i] - m;
        s += d*d;
    }
    s /= t;
    real log_likelihood = 0;
    for (int i=0; i<T; ++i) {
        log_likelihood += log_pdf(x[i]);
    }
    return log_likelihood;
}


/// Set Maximum Likelihood parameters
///
/// x is normal distributed, but we care about the log-normal density.
///
/// Returns the log-likelihood of the parameters
real NormalDistribution::setMaximumLikelihoodParametersLogNormal (std::vector<real>& x)
{
    int T = x.size();
    real t = (real) T;
    m = Sum(x) / t;
    s = 0;
    for (int i=0; i<T; ++i) {
        real d = x[i] - m;
        s += d*d;
    }
    s /= t;
    real log_likelihood = 0;
    for (int i=0; i<T; ++i) {
        log_likelihood += log_pdf(x[i]) - x[i];
    }
    return log_likelihood;
}


//------------------ Unknown mean -----------------------//
real NormalDistributionUnknownMean::generate()
{
    return prior.generate();
}

real NormalDistributionUnknownMean::generate() const
{
    return prior.generate();
}

/** Marginal pdf of a normal distribution with unknown mean.
    
  TODO Check that the marginal likelihood is correctly calculated
*/
real NormalDistributionUnknownMean::marginal_pdf(real x) const
{

    real mean = mu_n / tau_n;
    real sigma = 1.0 / tau_n + 1.0 / tau;
    real d = (mean - x )/sigma;
    return exp(-0.5 * d*d)/(sqrt(2.0 * M_PI) * sigma);
}

real NormalDistributionUnknownMean::pdf(real x) const
{
    return prior.pdf(x);
}

void NormalDistributionUnknownMean::calculatePosterior(real x)
{
    n++;
    mu_n += tau * x;
    tau_n += tau;
    prior.setMean(mu_n);
    prior.setSTD(1.0 / tau_n);
    observations.setMean(mu_n);
}

real NormalDistributionUnknownMean::getMean() const
{
    return mu_n / tau_n;
}

//------------------- Unknown mean and precision -------------------------//
NormalUnknownMeanPrecision::NormalUnknownMeanPrecision()
    : marginal(1), marginal_mean(1)
{
    mu_0 = 0;
    tau_0 = 1.0;
    alpha_0 = 1.0;
    beta_0 = 1.0;
    Reset();
}

NormalUnknownMeanPrecision::NormalUnknownMeanPrecision(real mu_0_, real tau_0_)
    : marginal(1), marginal_mean(1), mu_0(mu_0_), tau_0(tau_0_)
{
    alpha_0 = 1;
    beta_0 = 1;
    Reset();
}
void NormalUnknownMeanPrecision::Reset()
{
    n = 0;
    sum = 0.0;
    tau_n = tau_0;
    mu_n = mu_0; 
    alpha_n = alpha_0;
    beta_n = beta_0;
    M_2n = 0;
    bx_n = 0;


    Matrix InvT(1);
    InvT(0,0) = 1.0 / beta_0;
    Vector mean(1);
    mean(0) = mu_0;

    marginal_mean.setDegrees(alpha_0);
    marginal_mean.setLocation(mean);
    marginal_mean.setPrecision(tau_0 * alpha_0 * InvT);

    marginal.setDegrees(alpha_0);
    marginal.setLocation(mean);
    marginal.setPrecision((tau_0 / (tau_0 + 1.0)) * alpha_0 * InvT);

    //printf ("Univariate location: %f\n", mu_0);
}

NormalUnknownMeanPrecision::~NormalUnknownMeanPrecision()
{
}

real NormalUnknownMeanPrecision::generate()
{
    return marginal_mean.generate();
}

/// Generate from the posterior. Uses the ranlib implementation
real NormalUnknownMeanPrecision::generate() const
{
    return marginal.generate();
}

real NormalUnknownMeanPrecision::Observe(real x)
{
    real p = marginal_pdf(x);
    calculatePosterior(x);
    return p;
}

real NormalUnknownMeanPrecision::pdf(real x) const
{
    //real d = (x - mu_n);
    //real ln_p = 0.5 * log(tau_n / (2.0 * M_PI)) - 0.5 * tau_n * d * d;
    //return exp(ln_p);
    return marginal_mean.pdf(Vector(x));

}

real NormalUnknownMeanPrecision::log_pdf(real x) const
{
    return marginal_mean.log_pdf(Vector(x));
}


/** The marginal pdf of the observations.
    
    This is a student-t distribution.
*/
real NormalUnknownMeanPrecision::marginal_pdf(real x) const
{
    //printf ("# Univariate ");
    //marginal.Show();
    return marginal.pdf(Vector(x));
}


real NormalUnknownMeanPrecision::getMean() const
{
    return mu_n;
}


/** Calculate the posterior parameter distribution given a new real observation.

    This is the standard approach, but, instead of calculating
    \f[
    M_{2,n} = \sum_{i=1}^n (x_i - \bar{x}_n)^2
    \f]
    at each pass, we use Knuth's algorithm to calculate
    \f[
    M_{2,n} = M_{2,n-1} + (x_n - \bar{x}_{n})(x_n - \bar{x}_{n-1})
    \f]

 */
void NormalUnknownMeanPrecision::calculatePosterior(real x)
{
    n++;
    real bxn_prev = bx_n;
    real rn = (real) n;
    real irn = 1 / rn;
    bx_n = bx_n + (x - bx_n) * irn;
    M_2n = M_2n + (x - bx_n)*(x - bxn_prev);
    real delta_mean = (bx_n - mu_0);
    
    real gamma = (rn - 1) * irn;
    mu_n = mu_n * gamma + (1 - gamma) * x;
    tau_n += 1.0;

    alpha_n += 1.0;
    beta_n = beta_0 +  M_2n + delta_mean*delta_mean * tau_0 * rn / tau_n;
    


    Matrix InvT(1);
    InvT(0,0) = 1.0 / beta_n;
    Vector mean(1);
    mean(0) = mu_n;

    //printf ("# U: %f %f %f %f %f\n", tau_n, alpha_n, InvT(0,0), M_2n, mu_n);

    marginal_mean.setDegrees(alpha_n);
    marginal_mean.setLocation(mean);
    marginal_mean.setPrecision(tau_n * alpha_n  * InvT);

    marginal.setDegrees(alpha_n);
    marginal.setLocation(mean);
    marginal.setPrecision((tau_n / (tau_n + 1.0)) * alpha_n  * InvT );

    //marginal.Show();
}

/** Get the log-likelihood of some observations */
real NormalUnknownMeanPrecision::LogLikelihood(const std::vector<real>& x, int K) const
{
    real log_likelihood = LOG_ZERO;
    ExponentialDistribution accuracy_prior(1.0);
    for (int k=0; k<K; ++k) {
        real accuracy = accuracy_prior.generate();
        real variance = 1.0 / accuracy;
        NormalDistribution mean_prior(0.0, variance);
        real mean = mean_prior.generate();
        real log_p = 0.0;
        int n = x.size();
        NormalDistribution sample(mean, accuracy);
        for (int i=0; i<n; ++i) {
            log_p += sample.log_pdf(x[i]);
        }
        log_likelihood = logAdd(log_likelihood, log_p);
    }
    log_likelihood -= log(K);
    return log_likelihood;
}

/** Get the log-likelihood of some observations for Log-Normal observations */
real NormalUnknownMeanPrecision::LogLikelihoodLogNormal(const std::vector<real>& x, int K) const
{
    real log_likelihood = LOG_ZERO;
    ExponentialDistribution accuracy_prior(1.0);
    for (int k=0; k<K; ++k) {
        real accuracy = accuracy_prior.generate();
        real variance = 1.0 / accuracy;
        NormalDistribution mean_prior(0.0, variance);
        real mean = mean_prior.generate();
        real log_p = 0.0;
        int n = x.size();
        NormalDistribution sample(mean, accuracy);
        for (int i=0; i<n; ++i) {
            log_p += sample.log_pdf(x[i]) - x[i];
        }
        log_likelihood = logAdd(log_likelihood, log_p);
    }
    log_likelihood -= log(K);
    return log_likelihood;
}


void NormalUnknownMeanPrecision::Show() const
{
    printf("Normal-Gamma - ");
    
}

