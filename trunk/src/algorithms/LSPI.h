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
	Matrix A;
	Vector b;
	Vector w;
	RBFBasisSet* bfs;
	Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples;
	FixedContinuousPolicy policy;
public:	
	LSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, RBFBasisSet* bfs_, Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples_);
	~LSPI();
	
	Vector BasisFunction(Vector state, int action);
	void LSTDQ();
	void LSTDQ_Fast();
	void PolicyIteration();
	void Reset();
	real getValue(Vector state, int action);
	FixedContinuousPolicy& ReturnPolicy()
	{
		return policy;
	}
};

#endif
