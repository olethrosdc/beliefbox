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

#ifndef KNN_MODEL_H
#define KNN_MODEL_H

#include "Vector.h"
#include "KDTree.h"
#include <list>
#include <vector>

class TrajectorySample
{
public:
    Vector s;
    int a;
    real r;
    Vector s2;
    real V;
    TrajectorySample(Vector s_, int a_, real r_, Vector s2_)
        : s(s_), a(a_), r(r_), s2(s2_), V(0.0)
    {
    }
        
    
};

class KNNModel
{
protected:
    int n_actions;
    int n_dim;
    std::vector<KDTree<TrajectorySample>*> kd_tree;
    //RBFBasisSet basis;
    std::list<TrajectorySample> samples;
    real gamma;
    bool optimistic_values;
public:	
    KNNModel(int n_actions, int n_dim, real gamma_ = 0.9);
    ~KNNModel();
    void AddSample(TrajectorySample x);
    void GetExpectedTransition(Vector& x, int a, real& r, Vector& y, int K, real b);
    void GetExpectedTransitionDiff(Vector& x, int a, real& r, Vector& y, int K, real b);
    real GetExpectedActionValue(Vector& x, int a, int K, real b);
    real GetExpectedValue(Vector& x, int K, real b);
    int GetBestAction(Vector& x, int K, real b);
    void UpdateValue(TrajectorySample& start_sample, int K, real b);
    void ValueIteration(int K, real b);
};




#endif
