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
#include "MultinomialDistribution.h"

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
    assert(N==(uint) values.size());
    Pr.resize(N);

    // find the parents of each node in the graph
    for (unsigned int n=0; n<N; ++n) {
        HalfEdgeListIterator e = graph.getFirstParent(n);
        // how many values can our conditioned variable take?
        int permutations = 1;
        // how many permutations from the conditioning variables?
        for (int p=0; p<graph.n_parents(n); ++p, ++e) {
            permutations *= values.size(e->node);
        }

        Pr[n].Resize(permutations, values.size(n));
        printf ("Making CPT for %d of size %d x %d\n", n, permutations, values.size(n));
        for (int i=0; i<permutations; ++i) {
            for (int j=0; j<values.size(n); ++j) {
                Pr[n](i,j) = 1.0 / (real) values.size(n);
            }
        }
    }
    _calculate_depth();
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

void DiscreteBN::dotFile(const char* fname)
{
    FILE* fout = fopen(fname, "w");
    if (!fout) {
        fprintf(stderr, "Could not open file name %s for writing\n", fname);
        return;
    }

#if 0
   const char* colour[] ={
        "red",
        "green",
        "blue",
        "yellow",
        "magenta",
        "cyan",
        "black"};
#endif

    fprintf (fout, "digraph DBN {\n");
    fprintf (fout, "ranksep=2; rankdir=LR; \n");
    int n_variables = graph.n_nodes();
    for (int n=0; n<n_variables; n++) {
        HalfEdgeListIterator e = graph.getFirstParent(n);
        if (graph.n_parents(n)) {
            for (int p=0; p<graph.n_parents(n); ++p, ++e) {
                fprintf (fout,
                         "x%d -> x%d\n",
                         e->node, n);
            }
        } else {
            fprintf (fout,
                     "x%d\n",
                     n);
            
        }
    }
    fprintf (fout, "}\n");
    fclose(fout);
}

/** Generate a vector of samples
    
    Iterate through all depths, generating values for each one.
 */
Vector DiscreteBN::generate()
{
    int n_variables = graph.n_nodes();
    Vector x(n_variables);
    for (uint depth=0; depth<depth_list.size(); ++depth) {
        for (uint i=0; i<depth_list[depth].size(); ++i) {
            int node = depth_list[depth][i];
            HalfEdgeListIterator e = graph.getFirstParent(node);
            int index = 0;
            int permutations = 1;

            // how many permutations from the conditioning variables?
            //printf ("# Pr[X_%d | ", node);
            for (int p=0; p<graph.n_parents(node); ++p, ++e) {
                int val = (int) x[e->node];
                //printf ("X_%d=%d", e->node, val);
                index += permutations * val;
                permutations *= values.size(e->node);
            }
            //printf("] = ");
            //Vector probs = Pr[node].getRow(index);
            //for (int j=0; j<probs.Size(); ++j) {
            //printf ("%.2f ", probs[j]);
            //}
            //printf (" -- %d x %d\n", Pr[node].Rows(), Pr[node].Columns());
            MultinomialDistribution p(Pr[node].getRow(index));
            x[node] = p.generateInt();
        }
    }
    return x;
}

/** The probability of the vector of variables obtaining a particular value.

 */
real DiscreteBN::getLogProbability(std::vector<int>& x)
{
    assert((int) x.size() == graph.n_nodes());
    real log_p = 0;
    for (uint depth=0; depth<depth_list.size(); ++depth) {
        for (uint i=0; i<depth_list[depth].size(); ++i) {
            int node = depth_list[depth][i];
            HalfEdgeListIterator e = graph.getFirstParent(node);
            int index = 0;
            int permutations = 1;

            // how many permutations from the conditioning variables?
            for (int p=0; p<graph.n_parents(node); ++p, ++e) {
                int val = (int) x[e->node];
                index += permutations * val;
                permutations *= values.size(e->node);
            }

            log_p += log(Pr[node](index, x[node]));
        }
    }
    return log_p;
}



/** Obtain the full joint distribution

    For a network with N variables, return a matrix with N + 1
    columns, where the first N columns are the values of those
    variabels and the last column is the probability of the variables
    obtaining those values jointly.
 */
Matrix DiscreteBN::getJointDistribution()
{  
    int n_variables = graph.n_nodes();
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
    }
    
    return P;
}

void DiscreteBN::_calculate_depth_rec(std::vector<uint>& depth, int node, uint d)
{
    if (d >= depth[node]) { 
        return;
    }
    depth[node] = d;
    HalfEdgeListIterator c = graph.getFirstChild(node);
    for (int i=0; i<graph.n_children(node); ++i, ++c) {
        _calculate_depth_rec(depth, c->node, depth[node] + 1);
    }
}

void DiscreteBN::_calculate_depth()
{
    int n_variables = graph.n_nodes();

    std::vector<uint> depth(n_variables);
    for (int n=0; n<n_variables; n++) {
        depth[n] = n_variables;
    }
    for (int n=0; n<n_variables; n++) {
        HalfEdgeListIterator e = graph.getFirstParent(n);
        if (graph.n_parents(n) == 0) {
            depth[n] = 0;
            HalfEdgeListIterator c = graph.getFirstChild(n);
            for (int i=0; i<graph.n_children(n); ++i, ++c) {
                _calculate_depth_rec(depth, c->node, depth[n] + 1);
            }
        }
    }

    int max_depth = Max(depth);
    depth_list.resize(max_depth + 1);
    for (uint d=0; d<depth_list.size(); ++d) {
        //depth_list[d].resize(0);
        printf ("Depth %d: ", d);
        for (int n=0; n<n_variables; ++n) {
            if (depth[n] == d) {
                depth_list[d].push_back(n);
                printf ("%d ", n);
            }
        }
        printf("\n");
    }

}
