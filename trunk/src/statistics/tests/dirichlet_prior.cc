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

int main (void)
{
    int N = 4;
    DirichletDistribution dirichlet(N);
	Vector pre = dirichlet.GetParameters();
	Vector data(N);
    Vector theta(N);
    for (int i=0; i<N; ++i) {
        theta(i) = urandom();
    }
    theta /= theta.Sum();
    MultinomialDistribution P(theta);

    int c = 100;
    for (int t=0; t<10000; t++) {
        Vector x = P.generate();

        dirichlet.update(&x);

        Vector post = dirichlet.GetParameters();
        Vector gen = dirichlet.generate();
        c--;
        if (c == 0) {
            for (int i=0; i<N; i++) {
                printf ("%d %f %f %f %f\n",
                        i, pre[i]/pre.Sum(),
                        theta[i], post[i]/post.Sum(), gen[i]);
            }
            c = 100;
        }
        pre = post;
    }

    return 0;
}

#endif
