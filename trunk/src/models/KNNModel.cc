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
#include "Distribution.h"

KNNModel::KNNModel(int n_actions_, int n_dim_, real gamma_)
    : n_actions(n_actions_), n_dim(n_dim_), kd_tree(n_actions_), gamma(gamma_)
{
    for (int i=0; i<n_actions; ++i) {
        kd_tree[i] = new KDTree<TrajectorySample> (n_dim);
    }
    optimistic_values = true;
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
    Vector dummy(n_dim);
    RBF rbf(dummy, b);

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



/// Predict the next reward and a state
///
/// This is the simplest type of predictor.
/// It just predicts the next vector mean, no distribution is used.
void KNNModel::GetExpectedTransitionDiff(Vector& x, int a, real& r, Vector& y, int K, real b)
{

    Vector dummy(n_dim);
    RBF rbf(dummy, b);

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
        y += (sample->s2 - sample->s)* w;
        r += sample->r * w;
        sum += w;
    }
    y /= sum;
    y += x;
    r /= sum;
}

real KNNModel::GetExpectedActionValue(Vector& x, int a, int K, real b)
{
    Vector dummy(n_dim);
    RBF rbf(dummy, b);

    OrderedFixedList<KDNode> node_list = kd_tree[a]->FindKNearestNeighbours(x, K);
    
    std::list<std::pair<real, KDNode*> >::iterator it;
    
    real sum = 0;
    real Q = 0.0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
        KDNode* node = it->second;
        TrajectorySample* sample = kd_tree[a]->getObject(node);
        rbf.center = sample->s;
        real w =  rbf.Evaluate(x);
        Q += sample->V * w;
        sum += w;
    }
    if (optimistic_values)
    {
        real w = 0.01;
        Q += 1.0 *w;
        sum += w;
    }
    Q /= sum;
    return Q;
}

// Return the expected value of a state according to our model
real KNNModel::GetExpectedValue(Vector& x, int K, real b)
{


    Vector Q(n_actions);
    for (int a=0; a<n_actions; a++) {
        Q[a] = GetExpectedActionValue(x, a, K, b);

    }
    return Max(&Q);
}

// Return the expected value of a state according to our model
int KNNModel::GetBestAction(Vector& x, int K, real b)
{

    Vector Q(n_actions);
    Vector p(n_actions);
    for (int a=0; a<n_actions; a++) {
        Q[a] = GetExpectedActionValue(x, a, K, b);

    }
    SoftMax(Q.Size(), &Q[0], &p[0], 1.0);
    return DiscreteDistribution::generate(p);
}


/// Update the value of the current state
void KNNModel::UpdateValue(TrajectorySample& start_sample, int K, real b)
{
    Vector Q(n_actions);
    for (int a=0; a<n_actions; ++a) {
        RBF rbf(start_sample.s, b);
        OrderedFixedList<KDNode> node_list = kd_tree[a]->FindKNearestNeighbours(start_sample.s, K);
        
        std::list<std::pair<real, KDNode*> >::iterator it;
    
        real sum = 0;
        Q[a] = 0.0;
        for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
            KDNode* node = it->second;
            TrajectorySample* sample = kd_tree[a]->getObject(node);
            rbf.center = sample->s;
            real w =  rbf.Evaluate(start_sample.s);
            Q[a] += sample->r + gamma*GetExpectedValue(sample->s2, K, b);
            sum += w;
        }
        Q[a] /= sum;
    }
    start_sample.V = Max(&Q);
    
}

void KNNModel::ValueIteration(int K, real b)
{
    for (std::list<TrajectorySample>::iterator it = samples.begin();
         it != samples.end(); ++it) {
        UpdateValue(*it, K, b);
    }
    
}
