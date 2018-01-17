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

/// Constructor
KNNRegression::KNNRegression(int m, int n) : M(m), N(n), kd_tree(m) 
{
}

/// Add an element
void KNNRegression::AddElement(const PointPair& p)
{
    pairs.push_back(p);
    kd_tree.AddVectorObject(p.x, &pairs.back());
    //basis.AddCenter(p.x, 1.0);
}

/// Obtain a K-nearest neighbour estimate of E[y | x].
void KNNRegression::Evaluate(const Vector& x, Vector& y, int K) const
{

    RBF rbf(x, 10.0);

    //basis.Evaluate(x);
    assert(N == y.Size());
    for (int i=0; i<N; ++i) {
        y[i] = 0;
    }

    OrderedFixedList<KDNode> node_list = kd_tree.FindKNearestNeighbours(x, K);
    
    std::list<std::pair<real, KDNode*> >::iterator it;
    
    real sum = 0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
        KDNode* node = it->second;
        const PointPair* point_pair = kd_tree.getObject(node);
        real w = rbf.Evaluate(point_pair->x);
		//printf("R: "); point_pair->x.print(stdout);
		//printf("X: "); rbf.center.print(stdout);
		//printf ("%f ", w); (point_pair->y).print(stdout);
        y += point_pair->y * w;
        sum += w;
    }
	if (sum > 0) {
		y/=sum;
	}
	//printf("S: %f \n", sum);
    //y.print(stdout);
}



