/* -*- Mode: c++ -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "SparseGraph.h"
#include "ConnectivityMatrix.h"
#include "EasyClock.h"
#include <iostream>
#include <exception>
#include <stdexcept>

int main(int argc, char** argv)
{
  
    int N = 10;
    if (argc>=2) {
      N = atoi(argv[1]);
    }

    int K = 10;
    if (argc>=3) {
      K = atoi(argv[2]);
    }

    SparseGraph sparse(N, true);
    int parents = N/2;
    
    for (int k=0; k<K; ++k) {
        Edge e;
        e.src = rand()%parents;
        e.dst = rand()%parents + parents;
        printf ("%d -> %d\n", e.src, e.dst);
        e.w = 1.0;
        sparse.AddEdge(e);
    }
    
    ConnectivityMatrix matrix(sparse);
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            if (matrix.edge(i,j)^sparse.edge(i,j)) {
                std::cerr << "Inconsistent edges detected " << i << "->" << j << std::endl;
            }
        }
    }
    
    
    if (sparse.hasCycles()) {
        printf ("# ERR: SparseGraph should not have cycles\n");
    } else {
        printf ("# OK: no cycles detected in SparseGraph\n");
    }

    if (matrix.hasCycles()) {
        printf ("# ERR: ConnectivityMatrix should not have cycles\n");
    } else {
        printf ("# OK: no cycles detected in ConnectivityMatrix\n");
    }

    {
        {
            Edge e;
            e.src = N-1;
            e.dst = 0;
            e.w = 0;
            sparse.AddEdge(e);
        }
        for (int i=0; i<parents; i++) {
            Edge e;
            e.w = 1.0;
            e.dst = i+1;
            e.src = i;
            sparse.AddEdge(e);
            //printf ("%d -> %d\n", e.src, e.dst);
            e.dst = parents + i + 1;
            e.src = parents + i;
            if (e.dst < N) {
                sparse.AddEdge(e);
                //printf ("%d -> %d\n", e.src, e.dst);

            }
        }


        ConnectivityMatrix matrix(sparse);
        for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            if (matrix.edge(i,j)^sparse.edge(i,j)) {
                std::cerr << "Inconsistent edges detected " << i << "->" << j << std::endl;
            }
        }
        }
        if (sparse.hasCycles()) {
            printf ("# OK: SparseGraph should have cycles\n");
        } else {
            printf ("# ERR: no cycles detected in SparseGraph\n");
        }
        
        if (matrix.hasCycles()) {
            printf ("# OK: ConnectivityMatrix should have cycles\n");
        } else {
            printf ("# ERR: no cycles detected in ConnectivityMatrix\n");
        } 
    }
    return 0;
}

#endif
