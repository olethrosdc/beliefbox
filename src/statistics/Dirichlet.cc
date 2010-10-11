/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Dirichlet.h"
#include "ranlib.h"

/// Create a placeholder Dirichlet
DirichletDistribution::DirichletDistribution()
{
    n = 0;
}

/// Create a Dirichlet with uniform parameters
DirichletDistribution::DirichletDistribution(int n, real p)
{
    this->n = n;
    a.Resize(n);
    for (int i=0; i<n; ++i) {
        a[i] = p;
    }
}

/// Initialise parameters from a vector
DirichletDistribution::DirichletDistribution(Vector& x) : n(x.Size()), a(x)
{
    for (int i=0; i<n; ++i) {
        assert(a[i] >= 0);
    }
}

void DirichletDistribution::resize(int n, real p)
{
    this->n = n;
    a.Resize(n);
    for (int i=0; i<n; ++i) {
        a[i] = p;
    }
}

DirichletDistribution::~DirichletDistribution()
{
}

Vector DirichletDistribution::generate()
{
    Vector x(n);
    generate(x);
    return x;
}

void DirichletDistribution::generate(Vector& y)
{
        //Vector y(n);
    real sum = 0.0;
    for (int i=0; i<n; i++) {
        y[i] = gengam(1.0, a[i]);
        sum += y[i];
    }
    real invsum = 1.0 / sum;
     y *= invsum;
}


/** Dirichlet distribution
    Gets the parameters of a multinomial distribution as input.
*/
real DirichletDistribution::pdf(const Vector& x) const
{
    assert(x.Size() == n);

    real log_prod = 0.0;
    real sum = 0.0;
    for (int i=0; i<n; i++) {
        real xi = x[i];
        if (xi<0) {
            Swarning ("Got a negative value for x[%d]:%f\n", i, xi);
            return 0.0;
        }
        sum += xi;
        log_prod += xi * a[i];
    }
    if (fabs(sum-1.0f)>0.001) {
        Swarning ("Vector x not a distribution apparently: sum=%f.  Returning 0.\n", sum);
        return 0.0;
    }
    return exp(log_prod);
}


void DirichletDistribution::Observation(int i)
{
    a[i] += 1.0;
}

Vector DirichletDistribution::GetParameters() const
{
	return a;
}

Vector DirichletDistribution::GetMean() const
{
    Vector p = a;
    real S = a.Sum();
    if (S > 0) {
        p /= S;
    } else {
        real invs = 1.0 / (real) p.Size();
        for (int i=0; i<p.Size(); i++) {
            p[i] = invs;
        }
    }
    //printf ("sum: %f\n", p.Sum());
    assert(fabs(p.Sum()-1.0) < 0.0001);
	return p;
}
