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

#ifndef BINOMIAL_DISTRIBUTION_H
#define BINOMIAL_DISTRIBUTION_H

#include "Distribution.h"

class BinomialDistribution : public ParametricDistribution
{
public:
    real p;
    long t;
    real s;
    BinomialDistribution(){p=0.5f;t=1;}
    BinomialDistribution(real p, long t, real s=1.0);
    virtual ~BinomialDistribution() {}
    virtual void setMean(real mean) {p = mean;}
    virtual void setVariance(real var) {} ///< set variance in a magic way
    virtual real generate();
    virtual real pdf(real x) const;
    virtual real getMean() const {return s*t*p;}
    void setPeriod(long t);
};

#endif
