/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Matrix.h"
#include "KFold.h"
#include "MersenneTwister.h"

#include <vector>

/// Constructor
KFold::KFold(Matrix& data_, int K_)
	: data(data_),
	  K(K_),
	  T(data.Rows()),
	  N(data.Columns()),
	  assignment(T)
{

	MersenneTwisterRNG rng;
	rng.manualSeed(1801174128);
	std::vector<int> totals(K);
	for (int k=0; k<K; ++k) {
		totals[k] = 0;
	}
	for (int t=0; t<T; ++T) {
		int x = rng.discrete_uniform(K);
		assignment[t] = x;
		totals[x] ++;
	}
	bool flag = true;
	while (flag) {
		flag = false;
		for (int k=0; k<K; ++k) {
			if (totals[k] == 0) {
				flag = true;
				// choose a random, assigned part
				int t = rng.discrete_uniform(T);
				int swapped = assignment[t];
				assignment[t] = k;
				totals[k]++;
				totals[swapped]--;
			}
		}
	}
}

