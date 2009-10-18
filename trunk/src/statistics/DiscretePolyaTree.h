/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_POLYA_TREE_H
#define DISCRETE_POLYA_TREE_H

#include "DiscreteVariable.h"
#include "SparseGraph.h"
#include <vector>
#include "Vector.h"
#include "Matrix.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/

/** 
    A node in a Polya tree
 */
class PolyaNode
{
public:
    const int depth;
    const int n_values;
    std::vector<class PolyaNode*> nodes;
    std::vector<real> w;
    PolyaNode(int depth_, DiscreteVector& values_)
        : depth(depth_),
          n_values(values_.size(depth)),
          w(n_values)
    {
        for (int i=0; i<n_values; ++i) {
            w[i] = 0.5;
        }
        if (depth + 1 < values_.size()) {
            nodes.resize(n_values);
            for (int i=0; i<n_values; ++i) {
                nodes[i] = new PolyaNode(depth + 1, values_);
            }
        }

    }
    ~PolyaNode()
    {            
        for (uint i=0; i<nodes.size(); ++i) {
            delete nodes[i];
        }
        
    }
    void Observe(std::vector<int>& x)
    {
        int i = x[depth];
        w[i] += 1.0;
        if (nodes.size()) {
            nodes[i]->Observe(x);
        }
    }
    real getProbability(std::vector<int>& x)
    {
        int i = x[depth];
        real p = w[i] / Sum(w);
        if (nodes.size()) {
            p *= nodes[i]->getProbability(x);
        }
        return p;
    }

};


class DiscretePolyaTree
{
protected:
    DiscreteVector values;
    PolyaNode* root;
public:
    DiscretePolyaTree(DiscreteVector values_);
    ~DiscretePolyaTree();
    void Observe(std::vector<int>& x);
    void Observe(Vector& x);
    Matrix getJointDistribution();
    real getProbability(std::vector<int>& x);
};

#endif
