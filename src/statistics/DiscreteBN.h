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

#ifndef DISCRETE_BN_H
#define DISCRETE_BN_H

#include "DiscreteVariable.h"
#include "SparseGraph.h"
#include <vector>
#include "Vector.h"
#include "Matrix.h"

/** A Bayesian Network.

We have some discrete variables, \f$X\f$.  Each
variable has some parents.
*/
class DiscreteBN
{
protected:
    int n_variables;
    SparseGraph& graph; ///< The garph of dependencies between variables
    //std::vector<Matrix> Nc; ///< A matrix of outcome counts, one for each variable
    std::vector<Matrix> Pr; ///< A matrix of probabilities, one for each variable
    DiscreteVector values; ///< A vector of possible values that each variable can take
    std::vector<std::vector<uint> > depth_list; ///< a vector opf nodes at each depth.

    void _calculate_depth();  ///< calculate the depth of all nodes
    void _calculate_depth_rec(std::vector<uint>& depth, int node, uint d); ///< recursionfor calculating the depth of all nodes
 public:
    DiscreteBN(DiscreteVector _values, SparseGraph& _graph);
    Matrix& getProbabilityMatrix(int n);
    void setProbabilityMatrix(int n, Matrix& P);
    void dotFile(const char* fname);
    void generate(std::vector<int>& x);
    std::vector<int> generate()
    {
        std::vector<int> x(n_variables);
        generate(x);
        return x;
    }
    real getLogProbability(std::vector<int>& x);
    real getProbability(std::vector<int>& x)
    {
        return exp(getLogProbability(x));
    }
    Matrix getJointDistribution();

};
#endif
