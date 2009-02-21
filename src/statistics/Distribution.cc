/* -*- Mode: C++; -*- */
// VER: $Id: Distribution.c,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Distribution.h"
#include "SmartAssert.h"
#include "ranlib.h"
#include "MersenneTwister.h"
#include "Random.h"


real UniformDistribution::generate()
{
    real x = min + urandom()*range;
    //printf ("Uniform (%f,%f): %f\n", min, range, x);
    return x;
}

real UniformDistribution::pdf(real x)
{
    if ((x >= min)&&(x <= min + range))
        return 1.0f/range;
    return 0.0f;
}

SingularDistribution::SingularDistribution(real m)
{
    this->m = m;
}

void SingularDistribution::setMean(real mean)
{
    m = mean;
}

real SingularDistribution::getMean()
{
	return m;
}

real SingularDistribution::getVariance()
{
    return 0.0;
}

real SingularDistribution::generate()
{
    return m;
}
real SingularDistribution::pdf(real x)
{
        if (x==m) {
            return 1.0f;
        } else {
            return 0.0f;
        }
}


BernoulliDistribution::BernoulliDistribution(real p)
{
    assert(p>=0.0f && p<=1.0f);
    this->p = p;
}


real BernoulliDistribution::generate()
{
    if (urandom()<p) {
        return 1.0f;
    }
    return 0.0f;
}

real BernoulliDistribution::pdf(real x)
{
    if (x==0.0f) {
        return 1.0f - p;
    } else if (x==1.0f) {
        return p;
    }
    return 0.0f;
}



real LaplacianDistribution::generate()
{
    real x = urandom(-1.0, 1.0);
    real absx = fabs (x);
    real sgnx;
    if (x>0.0) {
        sgnx = 1.0;
    } else {
        sgnx = -1.0;
    }
    
    return m + sgnx * log(1.0 - absx) / l;

}

real LaplacianDistribution::pdf(real x)
{
    return 0.5*l * exp (-l*fabs(x-m));
}

real ExponentialDistribution::generate()
{
    real x = urandom();
    return - log (1.0 - x) / l;
}

real ExponentialDistribution::pdf(real x)
{
    real d = x - m;
    if (d>0.0) {
        return l * exp (-l*d);
    }
    return 0.0;
}


DiscreteDistribution::DiscreteDistribution() {
    p = NULL; n_outcomes=0;
}

DiscreteDistribution::DiscreteDistribution(int N) {
    p = NULL;
    n_outcomes = 0;
    SMART_ASSERT (N>0)(N);
    p = (real*) malloc (sizeof(real) * N);
    n_outcomes = N;
    real invN = 1.0/((real) N);
    for (int i=0; i<N; i++) {
        p[i] = invN;
    }
}

Distribution* DiscreteDistribution::clone()
{
    DiscreteDistribution* d = new DiscreteDistribution(n_outcomes);
    for (int i=0; i<n_outcomes; ++i) {
        d->p[i] = p[i];
    }
    return d;
}


DiscreteDistribution::~DiscreteDistribution() {
    if (p) {
        free (p);
    }
}

real DiscreteDistribution::generate()
{
    real d=urandom();
    real sum = 0.0;
    for (int i=0; i<n_outcomes; i++) {
        sum += p[i];
        if (d < sum) {
            return (real) i;
        }
    }
    SMART_ASSERT (0)(sum);
    return 0.0;
}

int DiscreteDistribution::generate(std::vector<real> x) const
{
    real d=urandom();
    real sum = 0.0;
    int n = x.size();
    for (int i=0; i<n; i++) {
        sum += x[i];
        if (d < sum) {
            return  i;
        }
    }
    return rand()%n; 
}


real DiscreteDistribution::getMean()
{
    real sum = 0.0;
    for (int i=0; i<n_outcomes; i++) {
        sum += p[i] * (float) i;
    }
	return sum;
}

real DiscreteDistribution::pdf(real x)
{
    int i=(int) floor(x);
    if ((i>=0)&&(i<n_outcomes)) {
        return p[i];
    } 
    return 0.0;
}

