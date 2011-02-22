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
#include "Random.h"
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <stdexcept>

void SpeedTest();


int main()
{
    int N = 8;
    int n_errors = 0;
    Matrix I = Matrix::Unity(N,N);
    Matrix W(N,N);
    Matrix x(N,N);
    Vector u(N);
    Matrix D(N,1);

    for (int i=0; i<N; ++i) {
        u[i] = i + 1;
    }

    Matrix three_by_four(N + 1, N + 2);
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            W(i,j) = true_random(false);
            printf ("%f ", W(i,j));
        }
        D(i,0)=true_random(false);
        printf("\n");
    }
    
    
    Matrix X(N,N);
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            X(i,j) = true_random(false);
        }
    }

    printf("W:\n");
    W.print(stdout);

    printf("I:\n");
    I.print(stdout);
    
    printf("-----------\n");
    printf("# W==W*I: ");

    if (W!=I*W) {
        n_errors++;
        printf ("FAILED\n");
    } 
    if (W==I*W) {
        printf ("OK\n");
    }

    printf("# (W*W)' == (W*W')': ");
    if (W*Transpose(W) == Transpose(W*Transpose(W))) {
        printf("OK\n");
    } else {
        n_errors++; printf("FAILED\n");
    }

    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            X(i,j) = true_random(false);
        }
    }

    printf ("# I = 0 ");
    if (Matrix::Null(N,N)!=I) {
        printf ("OK"); 
    } else {
        n_errors++; printf ("FAILED"); 
    }

    bool caught;
    try {
        caught = false;
        printf("# Invalid matrix multiplication: ");
        three_by_four*I;
    }


    catch (std::domain_error) {
        caught = true;
        printf ("OK - Domain error caught\n");
    }
    
    if (!caught) {
        n_errors++;
        printf ("FAILED - Domain error not caught\n");
    }
    
 
    if (W*2.0 != W+W) {
        n_errors++;
        printf ("# matrix-scalar multiplication FAILED\n");
    } else {
        printf ("# matrix-scalar multiplication OK\n");
    }

    printf("u'*I: ");
    (Transpose(u)*I).print(stdout);
    printf ("# vector transposition and multiplication OK\n"); 

    printf("I*u: ");
    ((const Matrix&) I*(const Vector&) u).print(stdout);
    printf ("# OK\n"); 

    printf("V*z: ");
    Matrix V(2,4);
    Vector z(4);
    for (int i=0; i<V.Rows(); ++i) {
        for (int j=0; j<V.Columns(); ++j) {
            V(i,j) = true_random(false);
        }
    }
    for (int j=0; j<z.Size(); ++j) {
        z(j)= true_random(false);
    }
    ((const Matrix&) V*(const Vector&) z).print(stdout);

    try {
        caught = false;
        printf("u*I:\n");
        (u*I).print(stdout);
    }

    catch (std::domain_error) {
        caught = true;
        printf ("# OK: Domain error caught\n");
    }

    if (!caught) {
        n_errors++;
        printf ("# ERR - Domain error not caught\n");
    }

    {
        printf("u*D': ");
        Matrix X = (u*Transpose(D));
        printf ("# vector and matrix transpose OK\n"); 
    }

    try {
        caught = false;
        printf("u*D:\n");
        ((Matrix) (u*D)).print(stdout);
    }
    catch (std::domain_error) {
        caught = true;
        printf ("# OK: Domain error caught\n");
    }
    if (!caught) {
        n_errors++;
        printf ("# ERR - Domain error not caught\n");
    }



    printf("Q = W, W = 2W\n");
    Matrix Q = W;
    W = W*2;
    //W.print(stdout);

    printf ("# 2W - Q\n");
    //(W - Q*2).print(stdout);
    fflush(stdout);
    if (W-Q*2 != Matrix::Null(N,N)) {
        fprintf (stderr, "Difference not NULL\n");
        n_errors++;
    }

    if (n_errors) {
        printf ("%d tests FAILED\n", n_errors);
    } else {
        printf ("All tests OK\n");
    }
    
    printf("Testing Decomposition of matrix:\n");
    W = Q * Transpose(Q);
    W.print(stdout);
    printf("Cholesky:\n");
    W.Cholesky().print(stdout);
    
    printf("LU:\n");
    real det;
    std::vector<Matrix> LU= W.LUDecomposition(det);
    LU[0].print(stdout);
    LU[1].print(stdout);
    
    printf("Testing inverse:\n");
    W.print(stdout);
    printf("Inverse of W:\n");
    Matrix invW = W.Inverse();
    invW.print(stdout);
    printf("Inverse of Inverse of W:\n");
    Matrix reinvW = invW.Inverse();
    reinvW.print(stdout);
    Matrix ShouldBeUnity = W * invW;
    printf("Should be Unity:\n");
    ShouldBeUnity.print(stdout);

    printf ("All tests OK\n");
    printf ("Performance testing:\n");
    double start_time = GetCPU();
    SpeedTest();
    double end_time = GetCPU();

    printf ("Total time: %f\n", end_time - start_time);
    
    return n_errors;
}

void SpeedTest()
{
    for (int iter=0; iter<100; iter++) {

        int K = ceil(urandom() * 100);
        int M = ceil(urandom() * 100);
        int N = ceil(urandom() * 100);
        
        Matrix A(K, M);
        for (int i=0; i<K; ++i) {
            for (int j=0; j<M; ++j) {
                A(i,j) = urandom();
            }
        }
        Matrix B(M, N);
        for (int i=0; i<M; ++i) {
            for (int j=0; j<N; ++j) {
                B(i,j) = urandom();
            }
        }
        Matrix D(K, N);
        for (int i=0; i<K; ++i) {
            for (int j=0; j<N; ++j) {
                D(i,j) = urandom();
            }
        }

        for (int k=0; k<10; ++k) {
            D += A*B;
        }
        
    }
}




#endif

/*
          AMD 64 3200 | Intel Atom N270
REFERENCE 9.5 (10.5)  | 55 (56)
MULTIPLIC 7.5 (8.3)   | 62 (62)
*/
