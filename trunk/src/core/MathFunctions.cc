// -*- Mode: c++ -*-
/* VER: $Id: MathFunctions.c,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "MathFunctions.h"
#include <cmath>
#include <cassert>
#include <cstdio>




/// Softmax an array \f$p_i = e^{\beta Q_i}/\sum_j^n{e^{\beta Q_i}}\f$
void SoftMax (int n, real* Q, real* p, real beta)
{
    real sum = 0.0;

    // Necessary for avoiding infinities.
    int arg_max = ArgMax(n, Q);
    real max = Q[arg_max];

    for (int i=0; i<n; i++) {
        p[i] = exp (beta * (Q[i]-max));
        sum += p[i];
    }
    sum = 1.0/sum;
    for (int i=0; i<n; i++) {
        p[i] *= sum;
    }
}
    
/// Softmin an array \f$p_i = e^{-\beta Q_i}/\sum_j^n{e^{-\beta Q_i}}\f$
void SoftMin (int n, real* Q, real* p, real beta)
{
    real sum = 0.0;

    // Necessary for avoiding infinities.
    int arg_min = ArgMin(n, Q);
    real min = Q[arg_min];

    for (int i=0; i<n; i++) {
        p[i] = exp (-beta * (Q[i]-min));
        sum += p[i];
    }
    sum = 1.0/sum;
    for (int i=0; i<n; i++) {
        p[i] *= sum;
    }
}
    
    
/// Normalise a vector to a destination vector (low level)
///
/// src is the source vector.
/// dst is the destination vector.
/// n_elements is the number of elements. 
/// As pointers are raw, make sure n_elements is correct.
/// It is safe for src and dst to point at the same vector.
void Normalise (const real* src, real* dst, const int n_elements)
{
	printf ("!");
    real sum = 0;
    for (int i=0; i<n_elements; i++) {
        sum += src[i];
    }
    if (sum==0) {
        for (int i=0; i<n_elements; i++) {
            dst[i] = src[i];
        }
        return;
    }
    assert(sum>0);
	real isum = 1.0 / sum;
    for (int i=0; i<n_elements; i++) {
        dst[i] = src[i] * isum;
    }
}

/// Return \f$\sum_i^n |a_i-b_i|^2\f$
real SquareNorm (real* a, real* b, int n)
{
    real sum = 0;
    for (int i=0; i<n; i++) {
        register real d = (*a++) - (*b++);
        sum += d*d;
    }
    return sum;
}

/// Return \f$\left(\sum_i^n |a_i-b_i|^2\right)^{1/2}\f$
real EuclideanNorm (real* a, real* b, int n)
{
    register real sum = 0;
    for (int i=0; i<n; i++) {
        register real d = (*a++) - (*b++);
        sum += d*d;
    }
    return sqrt(sum);
}

/// Return \f$\left(\sum_i^n |a_i-b_i|^p\right)^{1/p}\f$
real LNorm (real* a, real* b, int n, real p)
{
    real sum = 0;
    for (int i=0; i<n; i++) {
        register real d = (*a++) - (*b++);
        sum += pow(d,p);
    }
    return pow(sum,1.0/p);
}

/// Return \f$\sum_i^n |a_i-b_i|\f$
real L1Norm (real* a, real* b, int n)
{
    real sum = 0;
    for (int i=0; i<n; i++) {
        register real d = (*a++) - (*b++);
        sum += fabs(d);
    }
    return sum;
}

/// Return \f$\sum_i^n a_i\f$
real Sum (real* a, int n)
{
    real sum = 0;
    for (register int i=0; i<n; i++) {
        sum += *a++;
    }
    return sum;
}


#ifdef USE_DOUBLE
#define MINUS_LOG_THRESHOLD -39.14
#else
#define MINUS_LOG_THRESHOLD -18.42
#endif

real logAdd(real x, real y)
{
    if (x < y) {
        real tmp = x;
        x = y;
        y = tmp;
    }

    real minusdif = y - x;
#ifdef DEBUG
    if (std::isnan(minusdif))
        fprintf (stderr, "LogAdd: minusdif (%f) y (%f) or x (%f) is nan",minusdif,y,x);
#endif
    if (minusdif < MINUS_LOG_THRESHOLD)
        return x;
    else
        return x + log1p(exp(minusdif));
}

real logSub(real x, real y)
{
    if (x < y)
        fprintf(stderr, "LogSub: x (%f) should be greater than y (%f)", x, y);

    real minusdif = y - x;
#ifdef DEBUG
    if (std::isnan(minusdif))
        fprintf(stderr, "LogSub: minusdif (%f) y (%f) or x (%f) is nan",minusdif,y,x);
#endif
    if (x == y)
        return LOG_ZERO;
    else if (minusdif < MINUS_LOG_THRESHOLD)
        return x;
    else
        return x + log1p(-exp(minusdif));
}


