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
#include "MultivariateNormal.h"
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

    NormalDistribution generating_distribution(mean, sqrt(variance));

	Vector lower_bound(2);
	Vector upper_bound(2);
	for (int i=0; i<2; ++i) {
		lower_bound(i) = -100;
		upper_bound(i) = 100;
	}
    //ContextTreeKDTree pdf(2, max_depth, lower_bound, upper_bound);
    Vector mu(1);
    Matrix S = Matrix::Unity(1,1);
    S(0,0) = 1.0;
    mu(0) = 1.23456789;
    real tau = 1.0;
    real alpha = 1.0;
    MultivariateNormal multivariate_normal(mu, S);

    //generating_distribution.Show();
    //multivariate_normal.Show();
#if 0
    for (real x=-10; x<10; x+=0.01) {
        Vector X(1);
        X(0) = x;
        printf ("%f %f %f #X\n",
                x,
                generating_distribution.pdf(x),
                multivariate_normal.pdf(X));
    }
#endif

    MultivariateNormalUnknownMeanPrecision normal_wishart(mu, tau, alpha, S);
    NormalUnknownMeanPrecision normal_gamma(mu(0), tau);

    int randomise = urandom()*100;
    for (int i=0; i<randomise; i++) {
        generating_distribution.generate();
    }
    real sum_log_ng = 0.0;
    real sum_log_nw = 0.0;
    for (int t=0; t<T; ++t) {
        Vector x(1);
		x(0) = generating_distribution.generate();
        real p_ng = normal_gamma.Observe(x(0));
        real p_nw = normal_wishart.Observe(x);
        printf ("%f %f %f\n", p_ng, p_nw, generating_distribution.pdf(x(0)));
        sum_log_ng += log(p_ng);
        sum_log_nw += log(p_nw);
    }
    printf ("# Log sum: %f %f\n", sum_log_ng, sum_log_nw);
#if 0
    for (int t=0; t<10; ++t) {
        printf ("%f %f\n",
                normal_gamma.generate(),
                normal_gamma.generate());
    }
#endif
    return 0;
}

#endif
