/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "GaussianProcess.h"

/// Create a new GP with observations in R^d
GaussianProcess::GaussianProcess(Matrix& Sigma_p_,
                                 real noise_variance_)
    : Sigma_p(Sigma_p_),
      noise_variance(noise_variance_),
      X2(Matrix::Null(Sigma_p.Rows(), Sigma_p.Columns()))
{
    Accuracy = Sigma_p.Inverse();
    A = Accuracy;
}
         
GaussianProcess::~GaussianProcess()
{
}

Vector GaussianProcess::generate()
{
    Serror ("Not implemented\n");
    exit(-1);
    //return Vector();
}    
         
real GaussianProcess::pdf(Vector& x, real y)
{
    return 0.0;
}

/// This implements weight space view of a  GP
void GaussianProcess::Observe(Vector& x, real y)
{
    // Update total covariance
    Matrix V;
    Product(&x, &x, &V);
    X2 += V;

    A = X2 / noise_variance + Accuracy;
    Matrix inv_A = A.Inverse();
    //mean = inv_A * X;
}


/// This implements functon space view of a GP
///
/// X has a number of columns equal to the amount of data.
void GaussianProcess::Observe(Matrix& X, Vector& y)
{
    int N = X.Columns();
    Matrix K(N, N);
    for (int i=0; i<N; ++i) {
        Vector S = K.getColumn(i);
        for (int j=0; j<N; ++j) {
            real delta = (S - K.getColumn(j)).SquareNorm();
            K(i,j) = exp(-0.5 * delta);
            if (i==j) {
                K(i,j) += noise_variance;
            }
        }
    }
    
    Matrix L = K.Cholesky();
    

}
