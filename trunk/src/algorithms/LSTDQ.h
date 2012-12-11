// -*- Mode: c++ -*-
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef LSTDQ_H
#define LSTDQ_H

#include "real.h"
#include "Vector.h"
#include "Matrix.h"
#include "BasisSet.h"
#include "ContinuousPolicy.h"
#include "RandomPolicy.h"
#include <vector>
#include "Demonstrations.h"

class LSTDQ
{
protected:
	real gamma;
	int n_dimension;
	int n_actions;
	int n_basis;
	int algorithm;
	Matrix A;
	Vector b;
	Vector w;
	RBFBasisSet& bfs;
	Demonstrations<Vector, int>& Samples;
	FixedContinuousPolicy policy;
public:	
	LSTDQ(real gamma_,
		 int n_dimension_,
		 int n_actions_,
		 RBFBasisSet& bfs_,
		 Demonstrations<Vector, int>& Samples_);

	LSTDQ(real gamma_,
		 int n_dimension_,
		 int n_actions_,
		 int algorithm_,
		 RBFBasisSet& bfs_,
		 Demonstrations<Vector, int>& Samples_);
	~LSTDQ();
	
	Vector BasisFunction(const Vector& state, int action) const;
	void Calculate();
	void Calculate_Opt();
	void Reset();
	real getValue(const Vector& state, int action) const;
	real getValue(const Vector& state) const
	{
		real V = getValue(state, 0);
		for (int i=1; i<n_actions; ++i) {
			real Vi = getValue(state, i);
			if (Vi > V) {
				V = Vi;
			}
		}
		return V;
	}
	
};

#endif
