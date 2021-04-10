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

#ifndef LSPI_H
#define LSPI_H

#include "real.h"
#include "Rollout.h"
#include "Vector.h"
#include "Matrix.h"
#include "BasisSet.h"
#include "ContinuousPolicy.h"
#include "RandomPolicy.h"
#include <vector>


class LSPI{
protected:
	real gamma;
	real Delta;
	int n_dimension;
	int n_actions;
	int n_basis;
	int max_iteration;
	int algorithm;
	Matrix A;
	Vector b;
	Vector w;
	BasisSet<Vector, int>& bfs;
	Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples;
	FixedContinuousPolicy policy;
public:	
	LSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, BasisSet<Vector, int>& bfs_, Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples_);
	LSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, int algorithm_, BasisSet<Vector, int>& bfs_, Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples_);
	~LSPI();
	
	Vector BasisFunction(const Vector& state, int action);
	void LSTDQ();
	void LSTDQ(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update = true);
	void LSTDQ_OPT();
	void PolicyIteration();
	void Reset();
	void Update();
	real getValue(const Vector& state, int action);
	FixedContinuousPolicy& ReturnPolicy()
	{
		return policy;
	}
};

#endif
