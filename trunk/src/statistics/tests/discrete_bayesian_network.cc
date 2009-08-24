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
    
    DiscreteBN dbn(variable_space, *graph);

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
    
    for (int t=0; t<n_samples; ++t) {
        Vector v = dbn.generate();
        for (int i=0; i<v.Size(); ++i) {
            printf("%d ",  (int) v[i]);
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
        int n_vars = (int) ceil(urandom(1,8));
        int max_values = 2;
        test_dbn(n_vars, max_values, 16);
    }
    return 0;
}

#endif
