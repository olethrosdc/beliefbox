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

/**
   \defgroup MathGroup Mathematical functions
*/
/*@{*/


/// Return the minimum number in an array
template<typename T>
T Min (int n, T* x)
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
int ArgMin (std::vector<T>& x)
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
int ArgMax (std::vector<T>& x)
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
real Min (std::vector<T>& x)
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
real Max (std::vector<T>& x)
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



void SoftMax (int n, real* Q, real* p, real beta);
void SoftMin (int n, real* Q, real* p, real beta);
void Normalise (real* src, real* dst, int n_elements);
real EuclideanNorm (real* a, real* b, int n);
real SquareNorm (real* a, real* b, int n);
real LNorm (real* a, real* b, int n, real p);
real Sum (real* a, int n);

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
