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
    Vector s; ///< starting state
    int a; ///< action taken
    real r; ///< reward received
    Vector s2; ///< next state
    real V; ///< value of state
    real dV; ///< value difference
    bool terminal; ///< if state is terminal state
    TrajectorySample(Vector s_, int a_, real r_, Vector s2_, bool terminal_=false)
        : s(s_), a(a_), r(r_), s2(s2_), V(0.0), dV(1.0), terminal(terminal_)
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
    real optimism;
    real r_max;
    int max_samples;
public:	
    KNNModel(int n_actions, int n_dim, real gamma_ = 0.9, bool optimistic = true, real optimism_=0.1, real r_max_=0.0);
    ~KNNModel();
    void AddSample(TrajectorySample x);
    void GetExpectedTransition(real alpha, Vector& x, int action, real& reward, Vector& y, int K, real b);
    real GetExpectedActionValue(Vector& x, int a, int K, real b);
    real GetExpectedValue(Vector& x, int K, real b);
    int GetBestAction(Vector& x, int K, real b);
    void UpdateValue(TrajectorySample& start_sample, real alpha, int K, real b);
    void ValueIteration(real alpha, int K, real b);
    void Show();
    void SetMaxSamples(int max_samples_)
    {
        max_samples = max_samples_;
    }
};




#endif
