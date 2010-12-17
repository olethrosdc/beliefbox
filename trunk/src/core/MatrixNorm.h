/* -*- Mode: c++ -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef MATRIX_NORM
#define MATRIX_NORM

#include "Matrix.h"

real MaxNorm(const Matrix& X);
real FrobeniusNorm(const Matrix& X);
real PNorm(const Matrix& X, const real p);

#endif
