/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2004 by Christos Dimitrakakis <dimitrak@idiap.ch>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include "Object.h"
#include "Vector.h"
#include "real.h"
#include <cmath>

/** The gamma function.
    
    \f$Gamma(x) = \int_0^\infty t^{x-1} e^{-t} \, dt\f$.
 */
inline float Gamma(float x)
{
    return tgammaf(x);
}

/** The gamma function, double version.
 */
inline double Gamma(double x)
{
    return tgamma(x);
}


/** The log gamma function.
 */
inline float logGamma(float x)
{
    return lgammaf(x);
}

/** The log gamma function, double version.
 */
inline double logGamma(double x)
{
    return lgamma(x);
}

/** The log two-variable Beta function */
real logBeta(real x, real y);
/** The two-variable Beta function */
real Beta(real x, real y);

/** The lector Beta function */
real logBeta(const Vector& x);
/** The vector Beta function */
real Beta(const Vector& x);

float betacf(float a, float b, float x);
float BetaInc(float x, float a, float b);

unsigned long binomial (int n, int k);


#endif
