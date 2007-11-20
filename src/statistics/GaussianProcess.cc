/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "GaussianProcess.h"

/// Create a new GP with observations in R^d
GaussianProcess::GaussianProcess(int d_) : N(0), d(d_)
{
}

/// Create a new GP with a set of initial kernels k
GaussianProcess::GaussianProcess(std::vector<Vector> k) : kernels(k)
{
    d = k[0].Size();
    N = k.size();
}

GaussianProcess::~GaussianProcess()
{
}

Vector GaussianProcess::generate()
{
    //return Vector();
}    

real GaussianProcess::pdf(Vector x)
{
    return 0.0;
}

void GaussianProcess::update(Vector* x)
{
}
