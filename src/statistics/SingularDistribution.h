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


#ifndef SINGULAR_DISTRIBUTION_H
#define SINGULAR_DISTRIBUTION_H

#include "Distribution.h"
#include "real.h"

/// \brief Singular distribution.
/// Special case.
class SingularDistribution : public ParametricDistribution
{
public:
    real m;
    SingularDistribution() {m =0.0f;}
    SingularDistribution(real m);
    virtual ~SingularDistribution() {}
    virtual void setVariance (real var) {}; ///< ignore variance
    virtual void setMean (real mean); ///< set the mean
    virtual real getMean () const;
    virtual real getVariance() const;
    virtual real generate() const;
    virtual real pdf(real x) const;
};

/// \brief Singular distribution.
/// Special case.
class UnknownSingularDistribution : public ConjugatePrior
{
public:
    const Distribution* prior;
    bool observed;
    SingularDistribution Q;
    UnknownSingularDistribution();
    UnknownSingularDistribution(const Distribution* prior_);
    virtual ~UnknownSingularDistribution();
    virtual void calculatePosterior(real x);
    virtual real Observe(real x);
    virtual real pdf(real x) const;
    virtual real marginal_pdf(real x) const;
    virtual real generate() const;
	virtual real generateMarginal() const;
	virtual real getMean() const;
};

#endif
