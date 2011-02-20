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
    //return Vector();
}    
         
real GaussianProcess::pdf(Vector& x, real y)
{
    return 0.0;
}

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
