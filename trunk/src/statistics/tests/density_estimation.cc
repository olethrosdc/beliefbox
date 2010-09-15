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
#include "ContextTreeKDTree.h"
#include "ContextTreeRealLine.h"
#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"


int main (int argc, char** argv)
{
    if (argc != 5) {
        Serror ("Usage: density_estimation T D Alpha Beta\n");
        exit(-1);
    }
    int T = atoi(argv[1]);
    if (T < 0) {
        Serror("T should be >= 0\n");
        exit(-1);
    }

    int max_depth = atoi(argv[2]);
    if (max_depth <= 0) {
        Serror("max_depth should be >= 0\n");
        exit(-1);
    }

    real Alpha = atof(argv[3]);
    if (Alpha < 0) {
        Serror("Alpha should be >= 0\n");
        exit(-1);
    }

    real Beta = atof(argv[4]);
    if (Beta < 0) {
        Serror("Beta should be >= 0\n");
        exit(-1);
    }

    //BetaDistribution distribution(Alpha,Beta);
    //BetaDistribution distribution2(2*Beta,Alpha);
    NormalDistribution distribution(Alpha,Beta);
    NormalDistribution distribution2(Beta,Alpha);

	Vector lower_bound(2);
	Vector upper_bound(2);
	for (int i=0; i<2; ++i) {
		lower_bound(i) = -10;
		upper_bound(i) = 10;
	}
    //ContextTreeKDTree pdf(2, max_depth, lower_bound, upper_bound);
    Vector mu(2);
    Matrix S = Matrix::Unity(2,2);
    real tau = 1.0;
    real alpha = 1.0;
    MultivariateNormalUnknownMeanPrecision pdf(mu, tau, alpha, S);

    int randomise = urandom()*100;
    for (int i=0; i<randomise; i++) {
        distribution.generate();
    }
    for (int t=0; t<T; ++t) {
        Vector x(2);
		x(0) = distribution.generate();
        x(1) = 0.5*x(0) + distribution2.generate();
        real p = pdf.Observe(x);
		//printf ("%f\n", p);
    }
#if 1
    for (real x=-10; x<10; x+=0.1) {
		for (real y=-10; y<10; y+=0.1) {
			Vector v(2);
			v[0] = x;
			v[1] = y;
			printf ("%f ", pdf.pdf(v));// distribution.pdf(x)*distribution2.pdf(y));
		}
		printf(" #Y\n");
	}
    pdf.Show();
#endif
    return 0;
}

#endif
