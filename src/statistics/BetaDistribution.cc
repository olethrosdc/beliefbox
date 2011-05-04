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
#include "BetaDistribution.h"
#include "SpecialFunctions.h"

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
	assert(alpha > 0 && beta >= 0 || alpha >= 0 && beta > 0);
    return genbet(alpha, beta);
}

/// Generate using ranlib
real BetaDistribution::generate() const
{
	assert(alpha > 0 && beta >= 0 || alpha >= 0 && beta > 0);
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
