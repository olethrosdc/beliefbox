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

#include "DiscreteDBN.h"
#include <stdexcept>

/** Discrete Dynamic Bayesian Network
    
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
DiscreteDBN::DiscreteDBN(DiscreteVector _values, SparseGraph& _graph, int prior) : graph(_graph), values(_values)
{
    if (graph.hasCycles()) {
        throw std::domain_error("A Bayesian ");
    }
    unsigned int N = graph.n_nodes();
    assert(N==(uint) values.size());
    n_parents.resize(N);
    Nc.resize(N);

    // find the parents of each node in the graph
    for (unsigned int n=0; n<N; ++n) {
        n_parents[n] = graph.n_parents(n);
        HalfEdgeListIterator e = graph.getFirstParent(n);
        // how many values can our conditioned variable take?
        int permutations = values.size(n);
        // how many permutations from the conditioning variables?
        for (int p=0; p<n_parents[n]; ++p, ++e) {
            permutations *= values.size(e->node);
        }
        Nc[n].resize(permutations);
        for (unsigned int i=0; i<Nc[n].size(); ++i) {
            Nc[n][i] = prior;
        }
    }
}

/** Just get an index so that we know which counter to increment
 */
int DiscreteDBN::getCountIndex(int n, std::vector<int> X)
{
    //n_parents[n] = graph.n_parents(n);
    assert((uint) values.size()==X.size());
    HalfEdgeListIterator e = graph.getFirstParent(n);
    
    int index = X[n];
    int N = values.size(n);
    for (int p=0; p<n_parents[n]; ++p, ++e) {
        index += N * X[e->node];
        N *= values.size(e->node);
    }
    return index;
}

/** Now we observe some variables.
    
    X is the vector of variables and we need to update the
    corresponding counts.
    
    We just search for nodes which have a parent.
        
*/
int DiscreteDBN::observe(std::vector<int> X)
{
    int n_updates = 0;
    // Probably we can optimise this by having a list of
    // nodes with parents
    // There should be no loops!
    for (unsigned int n=0; n<X.size(); ++n) {
        if (n_parents[n]) {
            Nc[n][getCountIndex(n, X)]++;
            n_updates++;
        }
    }

    return n_updates;
}

/** Now we need to infer probabilities of events
 */
