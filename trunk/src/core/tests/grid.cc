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
    int n_dimensions = 3;
    int K = 5;

    Vector L(n_dimensions);
    Vector U(n_dimensions);
    Vector M(n_dimensions);
    for (int i=0; i<n_dimensions; ++i) {
        L[i] = 0; //urandom();
        U[i] = 1.0; //L[i] + urandom();
        M[i] = (L[i] + U[i]) * 0.5;
    }

    Grid grid(L, U);

    int n_errors = 0;
    int n_points = 100;

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


    EvenGrid even_grid(L, U, K);
    printf ("# Making even grid with %d subdivisions for %d intervals\n",
            K, even_grid.getNIntervals());
    std::vector<int> counts(even_grid.getNIntervals());
    for (int i=0; i<even_grid.getNIntervals(); ++i) {
        counts[i] = 0;
    }
    n_points=1000000;
    real points_per_interval = (real) n_points / (real) even_grid.getNIntervals();
    for (int t=0; t<n_points; ++t) {
        Vector x(n_dimensions);
        for (int i=0; i<n_dimensions; ++i) {
            real alpha = urandom();
            x[i] = alpha*U[i] + (1 - alpha)*L[i];
        }
        counts[even_grid.getInterval(x)]++;
    }
    real err = 0;
    for (int i=0; i<even_grid.getNIntervals(); ++i) {
        err += fabs((real) counts[i] - points_per_interval);
    }
    real err_prob = 0.01;
    
    printf ("# average error: %f < %f\n", err / (real) n_points, sqrt(log(1/err_prob) / (2 * points_per_interval)));
    

    if (!n_errors) {
        printf("# OK - Test successful\n");
    } else {
        fprintf (stdout, "%d / %d errors found\n", n_errors, n_points);
    }

    for (int i=0; i<even_grid.getNIntervals(); ++i){
        even_grid.getCenter(i).print(stdout);
    }

    return n_errors;
    
}

#endif
