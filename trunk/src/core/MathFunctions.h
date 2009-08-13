/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "real.h"
#include <cmath>
#include <vector>
#include <limits>

/**
   \defgroup MathGroup Mathematical functions
*/
/*@{*/


/// Return the minimum number in an array
template<typename T>
T Min (const int n, T* x)
{
    T min = x[0];
    for (int i=1; i<n; i++) {
        if (min>x[i]) {
            min = x[i];
        }
    }
    return min;
}

/// Return the maximum number in an array
template<typename T>
T Max (int n, T* x)
{
    T max = x[0];
    for (int i=1; i<n; i++) {
        if (max<x[i]) {
            max = x[i];
        }
    }
    return max;
}

template<typename T>
int ArgMin (int n, T* x)
{
    T min = x[0];
    int arg_min = 0;
    for (int i=1; i<n; i++) {
        if (min>x[i]) {
            min = x[i];
            arg_min = i;
        }
    }
    return arg_min;
}
/// Return the index of the maximum number in an array
template<typename T>
int ArgMax (int n, T* x)
{
    T max = x[0];
    int arg_max = 0;
    for (int i=1; i<n; i++) {
        if (max<x[i]) {
            max = x[i];
            arg_max = i;
        }
    }
    return arg_max;
}

template<typename T>
int ArgMin (const std::vector<T>& x)
{
    int n = x.size();
    T min = x[0];
    int arg_min = 0;
    for (int i=1; i<n; i++) {
        if (min>x[i]) {
            min = x[i];
            arg_min = i;
        }
    }
    return arg_min;
}

/// Return the index of the maximum number in an array
template<typename T>
int ArgMax (const std::vector<T>& x)
{
    int n = x.size();
    T max = x[0];
    int arg_max = 0;
    for (int i=1; i<n; i++) {
        if (max<x[i]) {
            max = x[i];
            arg_max = i;
        }
    }
    return arg_max;
}

template<typename T>
real Min (const std::vector<T>& x)
{
    int n = x.size();
    T min = x[0];
    for (int i=1; i<n; i++) {
        if (min>x[i]) {
            min = x[i];
        }
    }
    return min;
}
/// Return the index of the maximum number in an array
template<typename T>
real Max (const std::vector<T>& x)
{
    int n = x.size();
    T max = x[0];
    for (int i=1; i<n; i++) {
        if (max<x[i]) {
            max = x[i];
        }
    }
    return max;
}

/// Return the index of the maximum number in an array
template<typename T>
real L1Norm (const std::vector<T>& x)
{
    int n = x.size();
    T sum = fabs(x[0]);
    for (int i=1; i<n; i++) {
        sum += fabs(x[1]);
    }
    return sum;
}

template<typename T>
real Sum (const std::vector<T>& x)
{
    int n = x.size();
    T sum = x[0];
    for (int i=1; i<n; i++) {
        sum += x[i];
    }
    return sum;
}

template<typename T>
real Mean (const std::vector<T>& x)
{
    if (x.size() == 0) {
        //return 0;
        return std::numeric_limits<real>::quiet_NaN();
    }
    return Sum(x) / (real) x.size();
}

template<typename T>
real Span (const std::vector<T>& x)
{
    if (x.size() == 0) {
        //return 0;
        return std::numeric_limits<real>::quiet_NaN();
    }
    return Max(x) - Min(x);
}




void SoftMax (int n, real* Q, real* p, real beta);
void SoftMin (int n, real* Q, real* p, real beta);
void Normalise (real* src, real* dst, int n_elements);
real EuclideanNorm (real* a, real* b, int n);
real SquareNorm (real* a, real* b, int n);
real LNorm (real* a, real* b, int n, real p);
real L1Norm (real* a, real* b, int n);
real Sum (real* a, int n);

real logAdd(real x, real y);
real logSub(real x, real y);

inline bool approx_eq(real x, real y, real acc=10e-6)
{
    return (fabs(x - y) <= acc);
}

/// Get the sign of a number.
template<class T>
inline const T sign(const T& x)
{
    if (x>0) {
        return 1;
    } else if (x<0) {
        return -1;
    } else {
        return 0;
    }
} 

/*@}*/

#endif
