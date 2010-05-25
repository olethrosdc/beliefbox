/* -*- Mode: c++;  -*- */
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

#include "Grid.h"
#include "Vector.h"
#include "Random.h"

int main (void)
{
    int n_dimensions = 2;
    Vector L(n_dimensions);
    Vector U(n_dimensions);
    Vector M(n_dimensions);
    for (int i=0; i<n_dimensions; ++i) {
        L[i] = urandom();
        U[i] = L[i] + urandom();
        M[i] = (L[i] + U[i]) * 0.5;
    }

    Grid grid(L, U);

    int n_errors = 0;
    int n_points = 1000;

    for (int t=0; t<n_points; ++t) {
        Vector x(n_dimensions);
        int d = 1;
        int y = 0;
        for (int i=0; i<n_dimensions; ++i) {
            real alpha = urandom();
            real Z = urandom();
            if (Z < 0.5) {
                x[i] = alpha * L[i] + (1 - alpha)*M[i];
            } else {
                x[i] = alpha * U[i] + (1 - alpha)*M[i];
                y += d;
            }
            d <<= 1;
        }
        if (y != grid.getInterval(x)) {
            n_errors++;
        }
    }
    
    if (!n_errors) {
        printf("Test successful\n");
    } else {
        fprintf (stdout, "%d / %d errors found\n", n_errors, n_points);
    }
    return n_errors;
    
}

#endif
