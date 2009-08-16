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

#include "KNNModel.h"

KNNModel::KNNModel(int n_actions_, int n_dim_)
    : n_actions(n_actions_), n_dim(n_dim_), kd_tree(n_actions_)
{
    for (int i=0; i<n_actions; ++i) {
        kd_tree[i] = new KDTree<TrajectorySample> (n_dim);
    }
}

KNNModel::~KNNModel()
{
    for (int i=0; i< n_actions; ++i) {
        delete kd_tree[i];
    }
}

KNNModel::AddElement(TrajectorySample x)
{
    samples.push_back(x);
    kd_tree[x.a].AddVectorObject(x.s, &samples.back());
}

void KNNRegression::Evaluate(Vector& x, Vector& y, real& r, int K)
{

    RBF rbf(x, 1.0);

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
        PointPair* point_pair = kd_tree.getObject(node);
        rbf.center = point_pair->x;
        real w =  rbf.Evaluate(x);
        y += point_pair->y * w;
        sum += w;
    }
    y/=sum;
    
}



