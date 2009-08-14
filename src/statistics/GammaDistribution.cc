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

GammaDistribution::GammaDistribution() : alpha(1.0), beta(1.0)
{
}


GammaDistribution::GammaDistribution(real alpha_, real beta_) : alpha(alpha_), beta(beta_)
{
}


/// Gamma has a support \f$[0, \infty)\f$, with pdf
/// \f$x^{\alpha-1} \frac{\beta^\alpha \exp(-\beta x)}{\Gamma(\alpha)}\f$,
/// with $k\geq0$
real GammaDistribution::pdf(real x)
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
    return genbet(alpha, beta);
}
