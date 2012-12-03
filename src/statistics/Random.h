/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include "real.h"
#include <cmath>
#include <iostream>
#include "Object.h"
#include "Vector.h"
#include "MathFunctions.h"


/**
   \defgroup StatisticsGroup Statistics and probability
*/
/*@{*/

void setRandomSeed(unsigned int seed);
unsigned long lrandom();
real urandom();
real urandom(real min, real max);
Vector urandom(const Vector& min, const Vector& max);
/// Give a true random number
/// When blocking is true, then that takes 2-3 s per number.
/// When blocking is false, then that takse 12 us per number
real true_random(bool blocking=true);

unsigned long true_random_bits(bool blocking=true);
/*@}*/

#endif
