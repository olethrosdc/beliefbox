/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MultinomialDistribution.h"
#include "Random.h"

MultinomialDistribution::MultinomialDistribution(Vector p_)
    : p(p_)
{
}

MultinomialDistribution::MultinomialDistribution()
{
}

MultinomialDistribution::MultinomialDistribution(int n)
{
	p.Resize(n);
	real prior = 1.0 / (real) n;
	for (int i=0; i<n; i++) {
		p[i] = prior;
	}
}


MultinomialDistribution::MultinomialDistribution(std::vector<real> p_)
{
    int n = p_.size();
	p.Resize(n);
	for (int i=0; i<n; i++) {
		p[i] = p_[i];
	}
}

MultinomialDistribution::~MultinomialDistribution()
{
}

void MultinomialDistribution::Resize(int n)
{
	p.Resize(n);
	real prior = 1.0 / (real) n;
	for (int i=0; i<n; i++) {
		p[i] = prior;
	}
}


void MultinomialDistribution::generate(Vector* x)
{
	*x = generate();
}
void MultinomialDistribution::generate(Vector& x)
{
	x = generate();
}

int MultinomialDistribution::generateInt()
{
    real d=urandom();
    real sum = 0.0;
	int n = p.Size();

    for (int i=0; i<n; i++) {
        sum += p[i];
        if (d <= sum) {
            return i;
        }
    }
    return rand()%n;
}
Vector MultinomialDistribution::generate()
{
    int z = generateInt();
    int n = p.Size();
	Vector x(n);

	for (int i=0; i<n; i++) {
		if (i==z) {
			x[i] = 1.0;
		} else {
			x[i] = 0.0;
		}
	}
	return x;
}

/// Gives the pdf at any vector on R^n even though it is normally non-zero
/// only on a simplex.
real MultinomialDistribution::pdf(Vector* x) const
{
	int n = p.Size();
	real log_density = LOG_ONE;
	for (int i=0; i<n; i++) {
		log_density += p[i] + (*x)[i];
	}
	return exp(log_density);
}
real MultinomialDistribution::pdf(Vector& x) const
{
    return pdf(&x);
}


