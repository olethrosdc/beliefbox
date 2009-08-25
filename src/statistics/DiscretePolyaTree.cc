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

#include "DiscretePolyaTree.h"

DiscretePolyaTree::DiscretePolyaTree(DiscreteVector values_)
    : values(values_)
{
    // set up the parameters recursively
    root = new PolyaNode(0, values);
}

DiscretePolyaTree::~DiscretePolyaTree()
{
    delete root;
}

void DiscretePolyaTree::Observe(std::vector<int>& x)
{
    root->Observe(x);
}

void DiscretePolyaTree::Observe(Vector& x)
{
    std::vector<int> y(x.Size());
    for (uint i=0; i<y.size(); ++i) {
        y[i] = (int) x[i];
    }
    root->Observe(y);
}

Matrix DiscretePolyaTree::getJointDistribution()
{
    int n_variables = values.size();
    std::vector<int> x(n_variables);
    int permutations = 1;
    for (int i=0; i<n_variables; ++i) {
        x[i] = 0;
        permutations *= values.size(i);
    }
    Matrix P(permutations, n_variables + 1);
    bool flag = false;
    int r = 0;
    while(!flag) {
        for (int i=0; i<n_variables; i++) {
            P(r, i) = x[i];

        }
        P(r, n_variables) = getProbability(x);
        flag = values.permute(x);
        r++;
    }
    
    return P;
}

real DiscretePolyaTree::getProbability(std::vector<int>& x)
{
    return root->getProbability(x);
}
