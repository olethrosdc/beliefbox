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

struct TimeStatistics
{
    real combined_time;
    real search_time;
};
TimeStatistics knn_regression_test(int n_points, int K)
{
    std::vector<Vector> X(n_points);
    std::vector<Vector> Y(n_points);

    TimeStatistics time_statistics;
    
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
        knn_regression.AddElement(PointPair(X[i], Y[i]));
    }
    real end_time = GetCPU();
    time_statistics.combined_time = end_time - start_time;
    start_time = end_time;
    for (int i=0; i<n_points; ++i) {
        knn_regression.Evaluate(X[i], Z, K);
    }
    end_time = GetCPU();
    time_statistics.search_time = end_time - start_time;
    return time_statistics;

}
int main (void)
{
    for (int n_points=100; n_points<=100000; n_points*=2) {
        for (int K=2; K<=16; ++K) {
            TimeStatistics time_statistics = knn_regression_test(n_points, K);
            printf ("%d %d %f %f\n", n_points, K, time_statistics.combined_time, time_statistics.search_time);
            fflush(stdout);
        }
    }
    return 0;
}
#endif
