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

#ifndef NORMAL_DISTRIBUTION_H
#define NORMAL_DISTRIBUTION_H

#include "Distribution.h"
#include "Matrix.h"
#include "Vector.h"


/// Gaussian probability distribution
class NormalDistribution : public ParametricDistribution {
 private:
    bool cache;
    real normal_x, normal_y, normal_rho;
 public:
    real m; ///< mean
    real s; ///< standard deviation
    NormalDistribution() {m=0.0; s=1.0; cache = false;}
    /// Normal dist. with given mean and std
    NormalDistribution(real mean, real std)
    {
        setMean (mean);
        setVariance (std*std);
    }
    virtual Distribution* clone ()
    {
        NormalDistribution* d  = new NormalDistribution;
        d->m = m;
        d->s = s;
        return d;
    }
    virtual ~NormalDistribution() {}
    virtual real generate();
    virtual real pdf(real x);
    virtual void setVariance (real var) 
    {s = sqrt(var);} 
    virtual void setMean (real mean)
    {m = mean;}
    virtual real getMean ()
    {
        return m;
    }
};



/// Multivariate Gaussian probability distribution
class MultivariateNormal : public VectorDistribution {
 private:
    Vector mean;
    Matrix accuracy;
 public:

    MultivariateNormal(); 
    MultivariateNormal(Vector mean_, Matrix std_);
    virtual ~MultivariateNormal() {}
    virtual void generate(Vector& x);
    virtual Vector generate();
    virtual real pdf(Vector& x);
};

#endif
