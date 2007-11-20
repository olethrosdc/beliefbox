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

#ifndef DiscreteDBN_H
#define DiscreteDBN_H

#include "DiscreteVariable.h"
#include "SparseGraph.h"
#include <vector>
#include "Vector.h"

/** A Dynamic Bayesian Network.

We have some discrete variables, \f$A(t), X(t), X(t+1)\f$.  Each
variable has some parents such that \f$P(X(t+1) | X(t), A(t)) =
\prod_i P(X_i(t+1) | \phi_i)\f$.We count the number of times each
transition has been observed.
*/
class DiscreteDBN
{
protected:
    SparseGraph& graph;
    std::vector<int> n_parents;
    std::vector<std::vector<int> > Nc;
    DiscreteVector values;
    int getCountIndex(int n, std::vector<int> X);
public:
    DiscreteDBN(DiscreteVector _values, SparseGraph& _graph, int prior = 0);
    int observe(std::vector<int> X);
};
#endif
