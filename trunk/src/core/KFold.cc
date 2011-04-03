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
	  n_folds(K_),
	  n_records(data.Rows()),
	  n_columns(data.Columns()),
	  assignment(n_records),
	  totals(n_folds)
{
	MersenneTwisterRNG rng;
	rng.manualSeed(1801174128); 
	for (int k=0; k<n_folds; ++k) {
		totals[k] = 0;
	}
	for (int t=0; t<n_records; ++t) {
		int x = rng.discrete_uniform(n_folds);
		assignment[t] = x;
		totals[x] ++;
	}
	bool flag = true;
	while (flag) {
		flag = false;
		for (int k=0; k<n_folds; ++k) {
			if (totals[k] == 0) {
				flag = true;
				// choose a random, assigned part
				int t = rng.discrete_uniform(n_records);
				int swapped = assignment[t];
				assignment[t] = k;
				totals[k]++;
				totals[swapped]--;
			}
		}
	}
	printf("# Making %d-fold with totals: ", n_folds);
	for (int k=0; k<n_folds; ++k) {
		printf ("%d ", totals[k]);
	}
	printf("\n");
}

Matrix KFold::getTrainFold(int n)
{
	assert(n>=0 && n<n_folds);

	int M = n_records - totals[n];
	Matrix fold(M, n_columns);
	int m = 0;
	for (int t=0; t<n_records; ++t) {
		if (assignment[t]!=n) {
			for (int i=0; i<n_columns; ++i) {
				fold(m, i) = data(t, i);
			}
			m++;
		}
	}
	return fold;
}


Matrix KFold::getTestFold(int n)
{
	assert(n>=0 && n<n_folds);

	int M = totals[n];
	Matrix fold(M, n_columns);
	int m = 0;
	for (int t=0; t<n_records; ++t) {
		if (assignment[t]==n) {
			for (int i=0; i<n_columns; ++i) {
				fold(m, i) = data(t, i);
			}
			m++;
		}
	}
	return fold;
}
