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

#ifdef MAKE_MAIN
#include "DiscreteBN.h"
#include "Dirichlet.h"
#include "SparseGraph.h"
#include "DiscretePolyaTree.h"
#include "Random.h"

bool test_dbn(int n_variables, int max_values, int n_samples)
{
    if (max_values < 2) {
        fprintf (stderr, "No reason to do this with less than 2 values per variable!\n");
        exit(-1);
    }
    SparseGraph* graph = NULL;
    while (!graph) {
        graph = new SparseGraph(n_variables, true);
        int n_edges = (int) floor(urandom(0, 2*n_variables));
        for (int i=0; i<n_edges; ++i) {
            int src = rand()%n_variables;
            int dst = rand()%n_variables;
            if (src != dst && !graph->edge(src, dst) && !graph->edge(dst, src)) {
                Edge edge(src, dst);
                graph->AddEdge(edge, false);
                printf ("Adding edge %d -> %d\n", src, dst);
            }
        }
        if (graph->hasCycles()) {
            printf ("Cycles - trying again\n");
            delete graph;
            graph = NULL;
        }
    }
    std::vector<int> variable_space(n_variables);
    for (int i=0; i<n_variables; i++) {
        variable_space[i] = (int) floor(urandom(2.0, max_values+1.0));
    }
    DiscreteVector variable_specification(variable_space);
    
    DiscreteBN dbn(variable_specification, *graph);

    for (int n=0; n<n_variables; n++) {
        Matrix& P = dbn.getProbabilityMatrix(n);
        for (int i=0; i<P.Rows(); ++i) {
            real sum = 0.0;
            for (int j=0; j<P.Columns(); ++j) {
                real x= urandom();
                P(i,j) = x;
                sum += x;
            }
            for (int j=0; j<P.Columns(); ++j) {
                P(i,j) = P(i,j) / sum;
            }
        }
    }

    dbn.dotFile("dbn.dot");

    DiscretePolyaTree polya_tree(variable_specification);

    std::vector<int> v(n_variables);
    for (int t=0; t<n_samples; ++t) {
        dbn.generate(v);
        polya_tree.Observe(v);
        //for (uint i=0; i<v.size(); ++i) {
        ////            printf("%d ",  v[i]);
        //}
        //printf("# DATA\n");
    }
    
    Matrix P = dbn.getJointDistribution();
    Matrix Q = polya_tree.getJointDistribution();
    assert(P.Rows()==Q.Rows());
    for (int i=0; i<P.Rows(); ++i) {
        for (int j=0; j<P.Columns(); ++j) {
            printf ("(%f %f) ", P(i, j), Q(i, j));
        }
        printf("\n");
    }
    delete graph;

    return true;
}


int main (void)
{
    int n_iter = 1;
    for (int i=0; i<n_iter; ++i) {
        int n_vars = 3;//(int) ceil(urandom(1,8));
        int max_values = 2;
        test_dbn(n_vars, max_values, 163840);
    }
    return 0;
}

#endif
