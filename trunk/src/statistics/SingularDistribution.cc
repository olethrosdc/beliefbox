/* -*- Mode: C++; -*- */
// copyright (c) 2004-2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SingularDistribution.h"

SingularDistribution::SingularDistribution(real m)
{
    this->m = m;
}

void SingularDistribution::setMean(real mean)
{
    m = mean;
}

real SingularDistribution::getMean() const
{
	return m;
}

real SingularDistribution::getVariance() const
{
    return 0.0;
}

real SingularDistribution::generate() const
{
    return m;
}
real SingularDistribution::pdf(real x) const
{
        if (x==m) {
            return 1.0f;
        } else {
            return 0.0f;
        }
}

UnknownSingularDistribution::UnknownSingularDistribution(const Distribution& prior_) 
    : prior(prior_),
      observed(false)
{}

void UnknownSingularDistribution::calculatePosterior(real x)
{
    if (observed) {
        return;
    } 
}

real UnknownSingularDistribution::Observe(real x)
{
    if (observed) {
        return Q.pdf(x);
    }
    real p = prior.pdf(x);
    Q.m = x;
    return p;
}
