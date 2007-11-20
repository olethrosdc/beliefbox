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

    int K = N;
    if (argc==3) {
      K = atoi(argv[2]);
    }

    SparseGraph sparse(N, true);
    for (int k=1; k<K; ++k) {
        Edge e;
        e.src = rand()%N;
        e.dst = rand()%N;
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
    
    std::vector<float> distances(N);
    double start_time = GetCPU();
    for (int i=0; i<N; ++i) {
      matrix.CalculateDistances(&distances[0], i);
    }
    std::cout << "CM time " << GetCPU() - start_time << std::endl;
    start_time = GetCPU();
    for (int i=0; i<N; ++i) {
      sparse.CalculateDistances(&distances[0], i);
    }
    std::cout << "SP time " << GetCPU() - start_time << std::endl;
    return 0;
}

#endif
