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

#include "NormalDistribution.h"
#include "Random.h"

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
/// Normal distribution pdf
real NormalDistribution::pdf(real x) const
{
    real d = (m-x)/s;
    return exp(-0.5 * d*d)/(sqrt(2.0 * M_PI) * s);
}
real NormalDistributionUnknownMean::generate()
{
    return prior.generate();
}

real NormalDistributionUnknownMean::generate() const
{
    return prior.generate();
}

real NormalDistributionUnknownMean::pdf(real x) const
{
    return prior.pdf(x);
}

void NormalDistributionUnknownMean::calculatePosterior(real x)
{
    //    sum += x;
    n++;
    mu_n += tau * x;
    tau_n += tau;
    prior.setMean(mu_n);
    prior.setSTD(1.0 / tau_n);
    observations.setMean(mu_n);
}

/// Multivariate Gaussian generation
void MultivariateNormal::generate(Vector& x)
{
	x = mean;
}

/// Multivariate Gaussian generation
Vector MultivariateNormal::generate()
{	
	return mean; /// NOTE: FIX ME
}

/// Multivariate Gaussian density
real MultivariateNormal::pdf(Vector& x) const
{
	assert (x.Size()==mean.Size());
	real n = (real) x.Size();
	real d = (Transpose(x - mean) * accuracy * (x - mean))(0,0);
	real determinant = 1.0f; // well, not really.
	real log_pdf = - 0.5 * d - 0.5*(log(2.0*M_PI)*n + log(determinant));
	return exp(log_pdf);
}

