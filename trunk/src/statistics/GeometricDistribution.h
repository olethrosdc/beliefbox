/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "Distribution.h"

/// Geometric probability distribution
class GeometricDistribution : public ParametricDistribution {
public:
    real p; ///< parameter
    GeometricDistribution() 
    {
        setParameter(0.5);
    }
    /// Create an exponential distribution with parameter \c lambda
    GeometricDistribution(real p_)
    {
        setParameter(p_);
    }

    virtual ~GeometricDistribution() {}
    virtual real generate() const;
    virtual real pdf(real x) const;
    virtual void setVariance (real var)
    {
        assert(var > 0);
        setParameter((sqrt(1.0 + 4.0 * var) - 1) / (2.0 * var));
    }
    virtual void setMean (real mean)
    {
        assert(mean >= 0);
        setParameter(1.0 / (1.0 + mean));
    }
    virtual real getMean () const
    {
        return ((1.0 - p) / p);
    }
    void setParameter(real p_)
    {
        p= p_;
        assert (p >= 0 && p <= 1);
        gen_scale = 1.0 / log(1 - p);
    }
protected:
    real gen_scale; ///< scale for generation
    
};
