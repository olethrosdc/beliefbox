/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BinomialDistribution.h"
#include "SpecialFunctions.h"
#include "ranlib.h"

BinomialDistribution::BinomialDistribution(real p, long t, real s)
{
    this->p = p;
    this->t = t;
    this->s = s;
}

real BinomialDistribution::generate()
{
    return s * ignbin(t, p);
}

real BinomialDistribution::pdf(real x) const
{
    real logp = log(p);
    real log1mp = log(1-p);
    return binomial(t, (uint) x) * exp(logp*x + log1mp*(t-x));
}
void BinomialDistribution::setPeriod(long t)
{
    this->t = t;
}
