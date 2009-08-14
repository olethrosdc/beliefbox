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

#include "KNNRegression.h"

KNNRegression::KNNRegression(int n) : N(n), kd_tree(n) 
{
}
void KNNRegression::AddElement(PointPair p)
{
    pairs.push_back(p);
    kd_tree.AddVectorObject(p.x, &pairs.back());
    basis.AddCenter(p.x, 1.0);
}

void KNNRegression::Evaluate(Vector x, Vector& y)
{
    std::vector<real> d(point_set.size());
    for (uint i=1; i<point_set.size(); ++i) {
        x.distance(L);
    }
}



