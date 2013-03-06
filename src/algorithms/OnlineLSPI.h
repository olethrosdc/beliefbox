// -*- Mode: c++ -*-
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef Online_LSPI_H
#define Online_LSPI_H

#include "real.h"
#include "Rollout.h"
#include "Vector.h"
#include "Matrix.h"
#include "BasisSet.h"
#include "ContinuousPolicy.h"
#include "RandomPolicy.h"
#include <vector>

class OnlineLSPI{
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
	RBFBasisSet* bfs;
	FixedContinuousPolicy policy;
public:	
	OnlineLSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, RBFBasisSet* bfs_);
	OnlineLSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, int algorithm_, RBFBasisSet* bfs_);
	~OnlineLSPI();
	
	Vector BasisFunction(const Vector& state, int action);
	void LSTD(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update = true); 
	void LSTDQ(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update = true);
	void LSTDQ_OPT(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update = true);
	void Reset();
	void Update();
	real getValue(const Vector& state, int action);
	FixedContinuousPolicy& ReturnPolicy()
	{
		return policy;
	}
};

#endif
