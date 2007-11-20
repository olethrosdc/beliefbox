/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "Matrix.h"
#include "Distribution.h"
#include "NormalDistribution.h"
#include "EasyClock.h"
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <stdexcept>

void compute_determinants_both_ways(Matrix A, bool show = true)
{
    real det;
    if (show) {
        A.print(stdout);
    }
    Matrix B = A;

    double cpu_time = GetCPU();
    real det_sign = A.gaussian_elimination_forward();
    real det_gaus_el = det_sign * A.compute_det_triangular();
    cpu_time = GetCPU() - cpu_time;
    printf ("# det(A):%f (%f ms)\n", det_gaus_el, 1000.0*cpu_time);

    if (show) {
        printf ("# G:\n"); A.print(stdout);
    }
    cpu_time = GetCPU();
    std::vector<Matrix> LU = B.LUDecomposition(det);
    real det_lu_dec = det * LU[0].compute_det_triangular() * LU[1].compute_det_triangular();
    cpu_time = GetCPU() - cpu_time;
    if (show) {
        printf ("# L:\n"); LU[0].print(stdout);
        printf ("# U:\n"); LU[1].print(stdout);
    }

    printf ("# det(A)=%f (%f ms)\n", det_lu_dec, 1000.0*cpu_time);
    real delta = fabs(det_lu_dec - det_gaus_el);
    printf ("# delta = %f : ", delta);
    if (delta < 0.000001) {
        printf ("OK\n");
    } else {
        printf ("ERR\n");
    }
}

int main(void)
{
    Matrix Q(2,2);
    Matrix I = Matrix::Unity(2, 2);
    Q(0,0)=0.0; Q(0,1)=1.0;
    Q(1,0)=1.0; Q(1,1)=0.0;
    
    Matrix W(2,2);
    W(0,0)=1.0; W(0,1)=2.0;
    W(1,0)=3.0; W(1,1)=4.0;
    

    printf ("I\n");
    compute_determinants_both_ways(I);
    printf ("--------\n");
    printf ("Q\n");
    compute_determinants_both_ways(Q);
    printf ("--------\n");
    printf ("W\n");
    compute_determinants_both_ways(W);
    printf ("--------\n");


    for (int N=2; N<=1024; N*=2) {
        Matrix X = Matrix(N,N);
        NormalDistribution gaussian;
        real N2 = (real) (N);
        printf ("# N=%d\n", N);
        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                X(i,j) = gaussian.generate() / N2;
                if(i==j) X(i,j)+=1.0;
            }
        }
       
        compute_determinants_both_ways(X, false);
        printf ("------------------------------\n\n");
    }
    
    return 0;
}

#endif
