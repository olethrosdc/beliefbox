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

#ifndef EXPONENTIAL_DISTRIBUTION_H
#define EXPONENTIAL_DISTRIBUTION_H

#include "real.h"
class ParametricDistribution;

/// Exponential probability distribution
class ExponentialDistribution : public ParametricDistribution {
public:
    real l; ///< lambda
    real m; ///< mean
    ExponentialDistribution() {m=0.0; l=1.0;}
    /// Create an exponential distribution with parameter \c lambda
    ExponentialDistribution(real lambda)
    {
        l = lambda;
		m = 0.0;
    }
    ExponentialDistribution(real mean, real var)
    {
        setMean(mean);
        setVariance(var);
    }
    virtual ~ExponentialDistribution() {}
    virtual real generate() const;
    virtual real pdf(real x) const;
    virtual real log_pdf(real x) const;
    virtual void setVariance (real var)
    {l = sqrt(1.0f / var);}
    virtual void setMean (real mean)
    {m = mean;}
    virtual real getMean () const
    {
        return m;
    }
};


#endif
