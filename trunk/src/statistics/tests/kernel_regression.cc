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
#include "KernelRegression.h"
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

void kernel_regression_evaluate(int T)
{
	int n_x = 2;
	int n_y = 1;
	
	KernelRegression kernel_regression(n_x, n_y, 0.1, 5);

	Vector X(n_x);
	Vector Y(n_y);
	real theta = 0.0;

	for (int t=0; t<T; t++) {
		theta += 0.01;
		X(0) = sin(theta);
		X(1) = sin(theta + 0.01);
		Y(0) = 1.0;
        kernel_regression.AddPoint(X, Y);
		X(0) = sin(theta);
		X(1) = sin(theta);
		Y(0) = 0.0;
        kernel_regression.AddPoint(X, Y);
	}

	kernel_regression.BootstrapBandwidth();
	theta = 0;
	for (real x = -1; x<1; x+=0.01) {
		for (real y = -1; y<1; y+=0.01) {
			X(0) = x;
			X(0) = y;
			Vector Z = kernel_regression.expected_value(X);
			printf ("%f ", Z(0));
		}
		printf ("# theta\n");
	}
}

TimeStatistics kernel_regression_test(int n_points, int K)
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


    KernelRegression kernel_regression(n_in, n_out, 1.0, K);
    Vector Z(n_out);
    real start_time = GetCPU();
    for (int i=0; i<n_points; ++i) {
        //kernel_regression.Evaluate(X[i], Z, K);
        kernel_regression.AddPoint(X[i], Y[i]);
    }
    real end_time = GetCPU();
    time_statistics.combined_time = end_time - start_time;
    start_time = end_time;
    for (int i=0; i<n_points; ++i) {
        Vector Z = kernel_regression.expected_value(X[i]);
    }
    end_time = GetCPU();
    time_statistics.search_time = end_time - start_time;
    return time_statistics;

}
int main (void)
{
	kernel_regression_evaluate(100);
	
    for (int n_points=10; n_points<=100; n_points*=2) {
        for (int K=0; K<=16; ++K) {
            TimeStatistics time_statistics = kernel_regression_test(n_points, K);
            printf ("%d %d %f %f\n", n_points, K, time_statistics.combined_time, time_statistics.search_time);
            fflush(stdout);
        }
    }
    return 0;
}
#endif
