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
    DirichletDistribution dirichlet(2);
	Vector pre = dirichlet.GetParameters();
	Vector data(2);

    real actual_probability = 0.7;
    int c = 100;
    for (int t=0; t<10000; t++) {
        if (urandom() < actual_probability) {
            data[0] = 1.0;
            data[1] = 0.0;
        } else {
            data[0] = 0.0;
            data[1] = 1.0;
        }            

        dirichlet.update(&data);

        Vector post = dirichlet.GetParameters();
        Vector gen = dirichlet.generate();
        c--;
        if (c == 0) {
            for (int i=0; i<2; i++) {
                printf ("%d %f %f %f %f\n", i, pre[i], data[i], post[i], gen[i]);
            }
            c = 100;
        }
        pre = post;
    }

    return 0;
}

#endif
