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

/** K-Nearest-Neighbour regression */
class KNNRegression
{
protected:
    int M; ///< Tree and conditioning variable dimension
    int N; ///< Dimension of the conditioned variable
    KDTree<PointPair> kd_tree; ///< The tree
    //CoverTree<PointPair> kd_tree;
    //RBFBasisSet basis;
    std::list<PointPair> pairs; ///< A list of pairs
public:	
    KNNRegression(int m, int n);
    void AddElement(const PointPair& p);
    void Evaluate(Vector&x, Vector& y, int K);
};


#endif
