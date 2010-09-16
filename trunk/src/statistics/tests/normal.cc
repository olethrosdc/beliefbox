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
    if (argc !=  4) {
        Serror ("Usage: normal n_sample mean variance\n");
        exit(-1);
    }
    int T = atoi(argv[1]);
    if (T < 0) {
        Serror("T should be >= 0\n");
        exit(-1);
    }

    real mean = atof(argv[2]);
    if (mean < 0) {
        Serror("mean should be >= 0\n");
        exit(-1);
    }

    real variance = atof(argv[3]);
    if (variance < 0) {
        Serror("variance should be >= 0\n");
        exit(-1);
    }
    //BetaDistribution distribution2(2*Beta,Alpha);
    NormalDistribution distribution(mean, variance);

	Vector lower_bound(2);
	Vector upper_bound(2);
	for (int i=0; i<2; ++i) {
		lower_bound(i) = -100;
		upper_bound(i) = 100;
	}
    //ContextTreeKDTree pdf(2, max_depth, lower_bound, upper_bound);
    Vector mu(1);
    Matrix S = Matrix::Unity(1,1);
    real tau = 1.0;
    real alpha = 1.0;
    MultivariateNormalUnknownMeanPrecision normal_wishart(mu, tau, alpha, S);
    NormalUnknownMeanPrecision normal_gamma(mu(0), tau);

    int randomise = urandom()*100;
    for (int i=0; i<randomise; i++) {
        distribution.generate();
    }
    for (int t=0; t<T; ++t) {
        Vector x(2);
		x(0) = distribution.generate();
        normal_gamma.Observe(x(0));
        normal_wishart.Observe(x);
    }
#if 1
    for (real x=-10; x<10; x+=0.1) {
        Vector v(2);
        v[0] = x;
        printf ("%f %f %f %f", x, normal_wishart.pdf(v), normal_gamma.pdf(x), distribution.pdf(x));// distribution.pdf(x)*distribution2.pdf(y));
		printf(" #Y\n");
	}
    normal_wishart.Show();
#endif
    return 0;
}

#endif
