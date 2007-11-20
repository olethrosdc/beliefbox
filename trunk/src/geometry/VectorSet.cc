/* -*- Mode: C++; -*- */
/* $Revision$ */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "VectorSet.h"

VectorSet::VectorSet(int n)
{
    this->n = n;
}

VectorSet::~VectorSet()
{
}

void VectorSet::Add(Vector x)
{
    assert (x.n == n);
    points.push_back(x);
}

bool ConvexHull::Add(Vector x)
{
    if (IsInConvexHull(x)) {
        return true;
    }   
    return false;
}

bool ConvexHull::IsInConvexHull(Vector x)
{
    return false;
}
