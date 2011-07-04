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

/** Gamma distribution.
 *
 * The gamma distribution is conjugate with the exponential, Poisson
 * and the normal distribution for known mean.
 *
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
    virtual real generate();
    virtual real generate() const;
};

/// Conjugate prior to gamma distribution with unknown shape and
/// inverse scale
class GammaDistributionUnknownShapeScale : public ConjugatePrior
{
  real S; ///< sum of observations
  real P; ///< product of observations

}
#endif
