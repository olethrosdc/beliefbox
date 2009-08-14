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
#include "KNNRegression.h"
#include "Random.h"
#include <vector>
#include "Vector.h"
#include "Matrix.h"
#include "EasyClock.h"

bool knn_regression_test(int n_points, int K)
{
    std::vector<Vector> X(n_points);
    std::vector<Vector> Y(n_points);


    
    int n_in = 8;
    int n_out = 2;
    Matrix A(n_in, n_out);

    for (int j=0; j<n_in; ++j) {
        for (int k=0; k<n_out; ++k) {
            A(j,k) = urandom() - 0.5;
        }
    }    
    for (int i=0; i<n_points; ++i) {
        X[i].Resize(n_in);
        Y[i].Resize(n_out);

        for (int j=0; j<n_in; ++j) {
            X[i][j] = urandom();
        }

        for (int j=0; j<n_out; ++j) {
            Y[i][j] = 0.0;
        }

        for (int j=0; j<n_in; ++j) {
            real x = X[i][j];
            for (int k=0; k<n_out; ++k) {
                Y[i][k] += x * A(j, k);
            }
        }
    }


    KNNRegression knn_regression(n_in, n_out);
    Vector Z(n_out);
    real start_time = GetCPU();
    for (int i=0; i<n_points; ++i) {
        knn_regression.Evaluate(X[i], Z, K);
        real err = EuclideanNorm(&Y[i], &Z);
        printf ("%f\n", err);
        knn_regression.AddElement(PointPair(X[i], Y[i]));
    }
    real end_time = GetCPU();
    printf ("# TOTAL TIME: %f s \n", end_time - start_time);
    //for (int i=0; i<n_points; ++i) {
    //knn_regression.Evaluate(X[i], Z, K);
    //        printf ("%f %f %f\n", X[i][0], Y[i][0], Z[0]);
    //}

    return true;

}
int main (void)
{
    knn_regression_test(1000000, 5);
    return 0;
}
#endif
