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
#include "BasisSet.h"

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

void KNNModel::AddSample(TrajectorySample x)
{
    samples.push_back(x);
    kd_tree[x.a]->AddVectorObject(x.s, &samples.back());
}

/// Predict the next reward and a state
///
/// This is the simplest type of predictor.
/// It just predicts the next vector mean, no distribution is used.
void KNNModel::GetExpectedTransition(Vector& x, int a, real& r, Vector& y, int K, real b)
{

    RBF rbf(x, b);

    //basis.Evaluate(x);
    assert(n_dim == y.Size());
    r = 0.0;
    for (int i=0; i<n_dim; ++i) {
        y[i] = 0;
    }

    OrderedFixedList<KDNode> node_list = kd_tree[a]->FindKNearestNeighbours(x, K);
    
    std::list<std::pair<real, KDNode*> >::iterator it;
    
    real sum = 0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
        KDNode* node = it->second;
        TrajectorySample* sample = kd_tree[a]->getObject(node);
        rbf.center = sample->s;
        real w =  rbf.Evaluate(x);
        y += sample->s2 * w;
        r += sample->r * w;
        sum += w;
    }
    y /= sum;
    r /= sum;
}



