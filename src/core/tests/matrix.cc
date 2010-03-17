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
    int n_errors = 0;
    Matrix I = Matrix::Unity(2,2);
    Matrix W(2,2);
    Matrix x(2,2);
    Vector u(2);
    Matrix D(2,1);
    u[0] = 1;
    u[1] = 2;
    Matrix three_by_four(3,4);
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            W(i,j) = true_random(false);
            printf ("%f ", W(i,j));
        }
        D(i,0)=true_random(false);
        printf("\n");
    }
    
    
    Matrix X(2,2);
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
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

    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            X(i,j) = true_random(false);
        }
    }

    printf ("# I = 0 ");
    if (Matrix::Null(2,2)!=I) {
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
    if (W-Q*2 != Matrix::Null(2,2)) {
        fprintf (stderr, "Difference not NULL\n");
        n_errors++;
    }

    if (n_errors) {
        printf ("%d tests FAILED\n", n_errors);
    } else {
        printf ("All tests OK\n");
    }
    double start_time = GetCPU();
    SpeedTest();
    double end_time = GetCPU();
    printf ("Total time: %f\n", end_time - start_time);
    
    return n_errors;
}

void SpeedTest()
{
    for (int iter=0; iter<1000; iter++) {

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
