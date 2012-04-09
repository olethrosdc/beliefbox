/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MULTINOMIAL_DISTRIUBION_H
#define MULTINOMIAL_DISTRIUBION_H

#include "Vector.h"
#include "Distribution.h"

/// Multinomial distribution
class MultinomialDistribution : public VectorDistribution
{
protected:
	Vector p;
public:
	MultinomialDistribution(const std::vector<real>& p_);
	MultinomialDistribution(const Vector& p_);
	MultinomialDistribution(int n);
	MultinomialDistribution();
	virtual ~MultinomialDistribution();
	virtual void Resize(int n);
	virtual void generate(Vector* x) const;
	virtual Vector generate() const;
    int generateInt() const;
    virtual Vector getMean()
    {
        return p;
    }
    inline real& Pr(int i)
    {
        return p[i];
    }
    /// To enhance the API somewhat
	virtual real pdf(Vector* x) const
    {
        return pdf(*x);
    }
    virtual real pdf(const Vector& x) const;
    virtual void generate(Vector& x) const;

};

Vector MultinomialDeviation(const Vector& p, const int j, const real c);


#endif
