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

#include "MultivariateNormal.h"

//----------------------- Multivariate -----------------------------------//

MultivariateNormal::MultivariateNormal(const int n_dim_)
    : n_dim(n_dim_),
      mean(n_dim),
      accuracy(n_dim, n_dim)
{
    mean.Clear();
    accuracy.Clear();
    for (int i=0; i<n_dim; ++i) {
        accuracy(i,i) = 1;
    }
    determinant = 1;
}

MultivariateNormal::MultivariateNormal(const Vector& mean_, const Matrix& accuracy_)
    :  n_dim(mean_.Size()), mean(mean_), accuracy(accuracy_)
{
    accuracy.LUDecomposition(determinant);
}

/// In-place multivariate Gaussian generation
void MultivariateNormal::generate(Vector& x) const
{
	x = generate();
}

/** Multivariate Gaussian generation.

    Uses the Cholesky decomposition.
 */
Vector MultivariateNormal::generate() const
{	
    Matrix Sigma = accuracy;
    Sigma = Sigma.Inverse();
    Matrix A = Sigma.Cholesky();
    Vector v(n_dim);
    NormalDistribution normal;
    for (int i=0; i<n_dim; ++i) {
        v(i) = normal.generate();
    }
    const Matrix& Ar = A;
    return mean + Ar * v;
}


/** Multivariate Gaussian density.

    For a gaussian with mean and precision \f$\mu, T\f$, the pdf is given by
    \f[
    f(x \mid \mu, T) = (2\pi)^{-k/2} |T|^{1/2} 
    \f]
 */
real MultivariateNormal::log_pdf(const Vector& x) const
{
	assert (x.Size()==mean.Size());
	real n = (real) x.Size();
    Vector diff = x - mean;
	real d = Mahalanobis2(diff, accuracy, diff);
    assert(d >= 0);
	real log_pdf = 0.5 * (log(determinant) - d - n * log(2*M_PI));
    return log_pdf;
}


void MultivariateNormal::Show()  const
{
    printf("MultivariateNormal - ");
    printf("Mean: ");
    mean.print(stdout);
    printf("Accuracy:\n");
    accuracy.print(stdout);
}
