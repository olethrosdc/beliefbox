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
#include "Dirichlet.h"
#include "DirichletFiniteOutcomes.h"

int main (void)
{
    int N = 16384;
    DirichletDistribution dirichlet(N);
    DirichletFiniteOutcomes finite_dirichlet(N);

	Vector pre = dirichlet.GetParameters();
	Vector data(N);
    Vector theta(N);

    for (int i=0; i<N; ++i) {
        if (i >= 0) {
            theta(i) = 1.0; // (1.0 + (real) i);
        } else {
            theta(i) = 0.0;
        }
    }

    theta /= theta.Sum();
    MultinomialDistribution P(theta);

    int interval = 1000;
    int c = interval;
    for (int t=0; t<10000; t++) {
        Vector x = P.generate();

        dirichlet.update(&x);
        finite_dirichlet.update(&x);

        //Vector post = dirichlet.GetParameters();
        //Vector gen = dirichlet.generate();

        Vector post = dirichlet.getMarginal();
        Vector gen = finite_dirichlet.getMarginal();

        //Vector post = finite_dirichlet.getMarginal();
        //Vector gen = finite_dirichlet.generate();
        c--;
        if (c == 0) {
            real err1 = 0;
            real err2 = 0;
            for (int i=0; i<N; i++) {
                err1 += fabs(theta(i) - post(i));
                err2 += fabs(theta(i) - gen(i));
#if 1
                printf ("%d %f %f %f\n",
                        i,
                        theta(i),
                        post(i),
                        gen(i));
#endif
            }
            printf ("%f %f\n", err1, err2);
            c = interval;
        }
        pre = post;
    }

    return 0;
}

#endif
