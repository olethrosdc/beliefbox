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

#include "DiscreteBN.h"
#include <stdexcept>

/** Discrete Bayesian Network
    
    Arguments:
    _values: the specification of the discrete variable vector
    _graph: the dependencies between variables
    prior: the prior count

    Note that there is an trade-off between the number of variables
    and the number of paramters.  For example, if you have variable
    \f$Y\f$ conditioned on \f$X\f$, such that \f$Y_i \in \{0,1\}\f$
   depends on \f$X_j \in \{0,1,2,3\}\f$, then if f\$P(Y_i=1 | X_j=0)
    = q\f$ and \f$1-q\f$ for all other values of \f$X_j\f$ then one
    could instead create an auxilliary variable \f$Z = I(X_j=0)'f$ and
    connect \f$Y\f$ to that.  So we added one graph node and removed
    three parameters.
*/
DiscreteBN::DiscreteBN(DiscreteVector _values, SparseGraph& _graph) : graph(_graph), values(_values)
{
    if (graph.hasCycles()) {
        throw std::domain_error("A Bayesian network can not have cycles");
    }
    unsigned int N = graph.n_nodes();
    assert(N==values.size());
    Pr.resize(N);

    // find the parents of each node in the graph
    for (unsigned int n=0; n<N; ++n) {
        HalfEdgeListIterator e = graph.getFirstParent(n);
        // how many values can our conditioned variable take?
        int permutations = 1.0;
        // how many permutations from the conditioning variables?
        for (int p=0; p<graph.n_parents(n); ++p, ++e) {
            permutations *= values.size(e->node);
        }
        Pr[n].Resize(permutations, values.size(n));
        for (int i=0; i<permutations; ++i) {
            for (int j=0; j<values.size(n); ++j) {
                Pr[n](i,j) = 1.0 / (real) values.size(n);
            }
        }
    }
}

Matrix& DiscreteBN::getProbabilityMatrix(int n)
{
    assert(n >= 0 && n < graph.n_nodes());
    return Pr[n];
    
}
void DiscreteBN::setProbabilityMatrix(int n, Matrix& P)
{
    assert(n >= 0 && n < graph.n_nodes());
    assert(P.Columns() == Pr[n].Columns() && P.Rows() == Pr[n].Rows());
    Pr[n] = P;
}
