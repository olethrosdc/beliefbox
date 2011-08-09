/* -*- Mode: C++; -*- */
// copyright (c) 2004-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ranlib.h"
#include "BetaDistribution.h"
#include "SpecialFunctions.h"
#include "ExponentialDistribution.h"

/// Calculate
/// \f[
///  \exp(\log x(\alpha - 1) + \log (1-x)(\beta - 1) - B(\alpha, \beta)
///  = \frac {x^{\alpha -1}(1-x)^{\beta -1}}{B(\alpha, \beta)}
/// \f]
real BetaDistribution::pdf(real x) const
{
    if (x<0.0 || x>1.0) {
        return 0.0;
    }
    if (alpha == 1 && beta == 1) {
        return 1.0;
    }
    
    real log_pdf = 0;
    if (x == 0 && alpha == 1 && beta > 0) {
        log_pdf = -logBeta(alpha, beta);
    } else if (x == 1 && beta == 1 && alpha > 0) {
        log_pdf = -logBeta(alpha, beta);
    } else {
        log_pdf = log(x)*(alpha - 1.0) + log(1-x)*(beta - 1.0)- logBeta(alpha, beta);    
    }
    return exp(log_pdf);
}

real BetaDistribution::log_pdf(real x) const
{
    if (x<0.0 || x>1.0) {
        return log(0.0);
    }
    if (alpha == 1 && beta == 1) {
        return 0.0;
    }
    return log(x)*(alpha - 1.0) + log(1-x)*(beta - 1.0)- logBeta(alpha, beta);
}



/// Standard posterior calculation
void BetaDistribution::calculatePosterior(real x)
{
	assert (x>=0 && x <= 1);
    alpha += x;
    beta += (1.0-x);
}


void BetaDistribution::setMean(real mean)
{
    fprintf(stderr,"Warning: cannot set mean for Beta distribution\n");
} 

void BetaDistribution::setVariance(real var)
{
    fprintf(stderr, "Warning: cannot set variance for Beta distribution\n");
}

real BetaDistribution::getMean() const
{
    return alpha/(alpha + beta);
}

real BetaDistribution::getVariance() 
{
    real a_b = alpha + beta;
    return (alpha/a_b)*(beta/a_b)/(a_b + 1);
}

/// Generate using ranlib
real BetaDistribution::generate() 
{
	assert(alpha > 0 && (beta >= 0 || alpha >= 0) && beta > 0);
    return genbet(alpha, beta);
}

/// Generate using ranlib
real BetaDistribution::generate() const
{
	assert(alpha > 0 && (beta >= 0 || alpha >= 0) && beta > 0);
    return genbet(alpha, beta);
}

/// Generate using ranlib
real BetaDistribution::generateMarginal() 
{
	if (urandom() < getMean()) {
		return 1.0;
	} else {
		return 0.0;
	}
}

real BetaDistribution::setMaximumLikelihoodParameters(const std::vector<real>& x,
                                    int n_iterations)
{
    // First set up the mean and variance by the method of moments
    real Z = 1.0f / (real) x.size();
    real mean = Sum(x) * Z;

    real d = 0;
    for (uint i=0; i != x.size(); ++i) {
        assert(x[i] >= 0.0f && x[i] <= 1.0f);
        d += (x[i] - mean);
    }
    real variance = d * Z;
    real S = mean * (1.0f - mean) / variance - 1.0f;
    real max_alpha = 1;// mean * S;
    real max_beta = 1;//(1.0f - mean) * S;
    real max_log_likelihood = Distribution::log_pdf(x);
    //printf ("%f %f %f %d # alpha, beta, LL, iter\n", max_alpha, max_beta, max_log_likelihood, n_iterations);
    ExponentialDistribution Exp;
    for (int k=0; k<n_iterations; ++k) {
        alpha = 1 + Exp.generate();
        beta = 1 + Exp.generate();
        real log_likelihood = Distribution::log_pdf(x);
        if (log_likelihood > max_log_likelihood) {
            max_alpha = alpha;
            max_beta = beta;
            max_log_likelihood = log_likelihood;
            //printf ("%f %f %f # alpha, beta\n", max_alpha, max_beta, max_log_likelihood);
        }
    }
    alpha = max_alpha;
    beta = max_beta;
    return max_log_likelihood;
}
