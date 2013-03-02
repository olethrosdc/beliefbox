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
#include "MultivariateNormalUnknownMeanPrecision.h"
#include "BetaDistribution.h"
#include "iWishart.h"


int main (int argc, char** argv)
{
	//Vector mean(2);
//	Matrix sigma(2,2);
//	sigma(0,0) = 0.25;
//	sigma(0,1) = 0.3;
//	sigma(1,0) = 0.3;
//	sigma(1,1) = 1;
//	
//	int n = 31*31;
//	Matrix X(n,2);
//	int count = 0;
//	for(real i=-3.0; i<=3.1; i = i+0.2) {
//		for(real j=-3.0; j<=3.1; j = j + 0.2) {
//			printf("count = %d\n",count);
//			X(count,0) = i;
//			X(count,1) = j;
//			count++;
//		}
//	}
//	printf("dafdasdf\n");
//
//	MultivariateNormal m_normal(mean, sigma.Inverse());
//	Matrix F(n,1);
//	for(int i = 0; i<n; ++i) {
//		F(i,0) = m_normal.pdf(X.getRow(i));
//	}
//	
//	printf("X\n");
//	X.print(stdout);
//	printf("Result\n");
//	F.print(stdout);
    if (argc !=  4) {
        Serror ("Usage: normal n_sample scaling n_dimension\n");
        exit(-1);
    }
    int T = atoi(argv[1]);
    if (T < 0) {
        Serror("T should be >= 0\n");
        exit(-1);
    }

    real variance = atof(argv[2]);
    if (variance < 0) {
        Serror("variance should be >= 0\n");
        exit(-1);
    }

    int n_dim = atoi(argv[3]);
    if (n_dim <= 0) {
        Serror("n_dim must be > 0\n");
        exit(-1);
    }


    Vector mean(n_dim);
    for (int i=0; i<n_dim; ++i) {
        mean(i) = urandom();
    }
    logmsg("mean vector: "); mean.print(stdout);

    Matrix R = Matrix::Unity(n_dim, n_dim) / sqrt(variance);
    MultivariateNormal generating_distribution(mean, R);


    Vector mu(n_dim);
    Matrix S = Matrix::Unity(n_dim, n_dim);
    real tau = 1.0;
    real alpha = 1.0;

    MultivariateNormal multivariate_normal(mu, S);
    const Matrix& Sref = S;
    MultivariateNormalUnknownMeanPrecision normal_wishart(mu, tau, n_dim + alpha, Sref);

    real sum_log_nw = 0.0;
    for (int t=0; t<T; ++t) {
        Vector x = generating_distribution.generate();
        real p_nw = normal_wishart.Observe(x);
        //printf ("%f %f\n", p_nw, generating_distribution.pdf(x));
        sum_log_nw += log(p_nw);
    }
    printf ("# Log sum: %f\n", sum_log_nw);
#if 1
    for (int t=0; t<10; ++t) {
        printf ("x_%d %f\n", t,
                normal_wishart.generate()(0));
    }
#endif

#if 1
    // show the wishart
    normal_wishart.Show();

    // show the real precision matrix
    logmsg("Real precision:\n"); R.print(stdout);

    // generate covariance matrices from this wishart!
    iWishart wishart = normal_wishart.getPrecision();
	
    wishart.Show();

    for (int i=0; i<10; ++i) {
        Matrix A = wishart.generate();
        Matrix B = A.Inverse();
        printf("## %d %f %f %f\n", i,
               wishart.pdf(A),
               wishart.pdf(B),
               wishart.pdf(R.Inverse()));
        A.print(stdout);

        B.print(stdout);
    }
#endif
    return 0;
}

#endif
