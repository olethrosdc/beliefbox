/* -*- Mode: c++ -*- */
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

#include "Tensor.h"

#include <iostream>
#include <list>
#include <cstdlib>
using namespace std;

int main(void)
{
    int size=4;
    vector<int> X(size);
    
    for (int i=0; i<size; ++i) {
        X[i] = 1 + rand()%16;
    }

    Tensor tensor(X);
   
    cout << "This test should succeed\n";

    vector<int> x(size);
    int n_iter = 1000;
    for (int iter=0; iter<n_iter; ++iter) {
        for (int i=0; i<size; ++i) {
            x[i] = rand()%X[i];
        }
        tensor.Y(x) = iter;
    }

    cout << "This test should fail\n";
    for (int iter=0; iter<n_iter; ++iter) {
        for (int i=0; i<size; ++i) {
            x[i] = 1 + rand()%X[i];
        }
        tensor.Y(x) = iter;
    }
}

#endif
