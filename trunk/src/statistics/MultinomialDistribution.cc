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


/** Get the pdf.
    
    Note that 
    \f[
    f(x \mid p) = \prod_{i=1}^n p_i^{x_i},
    \f]
    if there is a single observation that is non-zero.
    Howeer, in general
    \f[
    f(x \mid p) = \|x\|_1! \prod_{i=1}^n p_i^{x_i} / {x_i!},
    \f]
 */

real MultinomialDistribution::pdf(const Vector& x) const
{
	int n = p.Size();
	real log_density = LOG_ONE;
    real sum_x = 0;
	for (int i=0; i<n; i++) {
		log_density += x(i) * log(p(i));
        sum_x += x(i);
	}
    if (sum_x >= 2) {
        Swarning("Unexpected value\n");
        for (int i=0; i<n; i++) {
            for (int k=0; k<x(i); k++) {
                log_density -= k;
            }
        }
        for (int k=0; k<sum_x; ++k) {
            log_density += k;
        }
    }
	return exp(log_density);
}


