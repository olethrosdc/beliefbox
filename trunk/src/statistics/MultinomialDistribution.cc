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

/// Construct it from a vector
MultinomialDistribution::MultinomialDistribution(const Vector& p_)
    : p(p_)
{
}

/// Create an empty one
MultinomialDistribution::MultinomialDistribution()
{
}

/// Create a uniform distribution on n outcomes
MultinomialDistribution::MultinomialDistribution(int n)
{
	p.Resize(n);
	real prior = 1.0 / (real) n;
	for (int i=0; i<n; i++) {
		p[i] = prior;
	}
}

/// Construct it from a std vector
MultinomialDistribution::MultinomialDistribution(const std::vector<real>& p_)
{
    int n = p_.size();
	p.Resize(n);
	for (int i=0; i<n; i++) {
		p[i] = p_[i];
	}
}

/// Destructor
MultinomialDistribution::~MultinomialDistribution()
{
}

/// resize to n elements and make uniform
void MultinomialDistribution::Resize(int n)
{
	p.Resize(n);
	real prior = 1.0 / (real) n;
	for (int i=0; i<n; i++) {
		p[i] = prior;
	}
}

/// generate from it
void MultinomialDistribution::generate(Vector* x) const
{
	*x = generate();
}


/// generate from it
void MultinomialDistribution::generate(Vector& x) const
{
	x = generate();
}

/// generate an integer
int MultinomialDistribution::generateInt(const Vector& x)
{
    real d=urandom();
    real sum = 0.0;
	int n = x.Size();

    for (int i=0; i<n; i++) {
        sum += x[i];
        if (d <= sum) {
            return i;
        }
    }
    return rand()%n;
}

/// generate and return a vector
Vector MultinomialDistribution::generate() const
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

/** Multinomial Deviation.

    Given a vector \f$p \in S^n\f$, an index \f$j \in N_n\f$, and a scalar \f$c \in R \f$, return the solution of:
    \f[
    \max \{\sign(c) (q_j - p_j) \mid q \in S^n, \|p - q\|_1 \leq |c|\}
    \]
*/
Vector MultinomialDeviation(const Vector& p, const int j, const real c)
{
    int n = p.Size();
    assert (n > 1);
    assert (approx_eq(p.Sum(), 1.0));
    assert (Min(p) >= 0.0 && Max(p) <= 1.0);
    
    real d_j = 0.5 * c;
    
    if (c > 0.0) {
        d_j = std::min(1 - p(j), d_j);
    } else {
        d_j = std::max(-p(j), d_j);
    }
    real d_rest = -d_j / (real) (n-1);
    Vector q(n);
    q +=  d_rest;
    q(j) = d_j;
    q += p;
    real s = 0.0;
    for (int i=0; i<n; ++i) {
        if (q(i) < 0.0) {
            q(i) = 0.0;
        }
        s += q(i);
    } 

    if (s > 1.0) {
        q /= s;
    }

    return q;
}
