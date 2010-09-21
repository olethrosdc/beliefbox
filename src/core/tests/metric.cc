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

#ifdef MAKE_MAIN

#include "Vector.h"
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <stdexcept>

real euclidean_metric(Vector& x, Vector& y)
{
	return EuclideanNorm(&x, &y);
}

int main()
{
	Vector x(2);
	Vector y(2);
	x[0] = 0.1;
	x[1] = 0.2;
	y[0] = 1;
	y[1]= 2;
	
	
}

#endif
