/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "GaussianProcess.h"
#include "NormalDistribution.h"

int main (int argc, char** argv)
{
    int T = 1000;
	int n_dim = 1;
	real noise_variance = 1.0;

	Matrix Sigma = Matrix::Unity(n_dim, n_dim);
	GaussianProcess GP(Sigma,
					   noise_variance);
	NormalDistribution normal;

	Matrix InputData(T, n_dim);
	for (int t=0; t<T; ++t) {
		for (int n=0; n<n_dim; ++n) {
			InputData(t, n) = normal.generate();
		}
	}
	
	Vector OutputData(T);
	for (int t=0; t<T; ++t) {
		real y = 0;
		for (int n=0; n<n_dim; ++n) {
			y += InputData(t, n) * ((real) t) * sin(n);
		}
	}
	
	GP.Observe(InputData, OutputData);
    return 0;
}

#endif
