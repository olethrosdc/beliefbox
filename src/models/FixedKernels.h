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

#include "Model.h"
#include "Vector.h"
#include <vector>

class TriangularKernel : public Model<real, real>
{
protected:
	real m, w;
public:

	TriangularKernel(real m_, real w_);
	virtual ~TriangularKernel() {};
	virtual real Get(real x);
};

class ConeKernel : public Model<Vector, Vector>
{
protected:
	Vector m;
	real w;
public:

	ConeKernel(Vector m_, real w_);
	virtual ~ConeKernel() {};
	virtual real Get(real x);
};


class TriangularKernelSet : public Model<real, real>
{
protected:
	std::vector<TriangularKernel> k;
	int n;
public:
	TriangularKernelSet(int n, real width, real lo, real hi);
	virtual ~TriangularKernelSet();
};

