/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "CoverTree.h"
#include "Random.h"
#include "Vector.h"
#include <vector>

int main()
{
    CoverTree tree;

    int n_points = 10;
    int n_dimensions = 1;
    printf ("Testing with %d points and %d dimensions\n", n_points, n_dimensions);
    std::vector<Vector> X(n_points);

    for (int i=0; i<n_points; i++) {
        X[i].Resize(n_dimensions);
        for (int j=0; j<n_dimensions; j++) {
            X[i][j] = (real) i / (real) n_points;//urandom();
        }
    }
    

    for (int i=0; i<n_points; i++) {
        tree.Insert(X[i]);
    }    

    return 0;
}

#endif
