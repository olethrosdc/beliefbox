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

/// General distributions package. Currently statistical estimation
/// is only provided through gradient descent methods for minimisation.

inline float Gamma(float x)
{
    return tgammaf(x);
}

inline double Gamma(double x)
{
    return tgamma(x);
}


inline float logGamma(float x)
{
    return lgammaf(x);
}

inline double logGamma(double x)
{
    return lgamma(x);
}

real logBeta(real x, real y);
real Beta(real x, real y);

real logBeta(Vector& x);
real Beta(Vector& x);

float betacf(float a, float b, float x);
float BetaInc(float x, float a, float b);


#endif
