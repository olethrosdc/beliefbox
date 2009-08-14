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

#ifndef KNN_REGRESSION_H
#define KNN_REGRESSION_H

#include "KDTree.h"
#include "PointPair.h"
#include "BasisSet.h"

class KNNRegression
{
protected:
    int M;
    int N;
    KDTree<PointPair> kd_tree;
    //RBFBasisSet basis;
    std::list<PointPair> pairs;
public:	
    KNNRegression(int m, int n);
    void AddElement(PointPair p);
    void Evaluate(Vector&x, Vector& y, int K);
};


#endif
