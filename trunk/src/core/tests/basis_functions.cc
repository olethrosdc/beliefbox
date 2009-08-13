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

#include "BasisSet.h"
#include "Random.h"

#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv)
{
    
    RBFBasisSet basis;
    int n_points = 8;
    real b = 10.0;
    std::vector<Vector> X(n_points);
    for (int i=0; i<n_points; ++i) {
        real x = (real) i / (real) n_points;
        X[i].Resize(1);
        X[i][0] = x;
        basis.AddCenter(X[i], b);
    }
    real dx = 0.01;
    for (real x = -1; x < 2.0; x += dx) {
        Vector V(1);
        V[0] = x;
        basis.Evaluate(V);
        printf("%f ", x);
        for (int i=0; i < basis.size(); i++) {
            printf ("%f ", basis.F(i));
        }
        printf("\n");
    }
}

#endif
