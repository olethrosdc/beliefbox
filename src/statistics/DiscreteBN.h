/* -*- Mode: C++; -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DiscreteBN_H
#define DiscreteBN_H

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
    SparseGraph& graph; ///< The garph of dependencies between variables
    //std::vector<Matrix> Nc; ///< A matrix of outcome counts, one for each variable
    std::vector<Matrix> Pr; ///< A matrix of probabilities, one for each variable
    DiscreteVector values; ///< A vector of possible values that each variable can take
 public:
    DiscreteBN(DiscreteVector _values, SparseGraph& _graph);
    Matrix& getProbabilityMatrix(int n);
    void setProbabilityMatrix(int n, Matrix& P);
};
#endif
