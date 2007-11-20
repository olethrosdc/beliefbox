// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "FixedKernels.h"
#include <cmath>

TriangularKernel::TriangularKernel(real m_, real w_) : m(m_), w(w_)
{
}


real TriangularKernel::Get(real x)
{
	real d = w - fabs(m-x);

	if (d < 0) {
		return 0.0;
	}
	return d;
}

TriangularKernelSet::TriangularKernelSet(int n_, real width, real lo, real hi) : n(n_)
{
	real d = (hi - lo) / (real) (n + 1);
	real x = lo;
	for (int i=0; i<n; ++i) {
		x+= d;
		k.push_back(TriangularKernel(x, width));
	}
}
