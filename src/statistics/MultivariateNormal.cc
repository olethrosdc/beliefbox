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

/// Multivariate Gaussian generation
void MultivariateNormal::generate(Vector& x) const
{
	x = mean;
}

/// Multivariate Gaussian generation
Vector MultivariateNormal::generate() const
{	
	return mean; /// NOTE: FIX ME
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
//----------------- Multivariate Unknown mean and precision -----------------------//
MultivariateNormalUnknownMeanPrecision::MultivariateNormalUnknownMeanPrecision()
    : n_dim(1), marginal_mean(1), marginal(1), bx_n(n_dim), M_2n(n_dim, n_dim)
{
    mu_0 = 0.0;
    tau_0 = 1.0;
    alpha_0 = 1;
    T_0 = Matrix::Unity(1,1);
    Reset();
}

MultivariateNormalUnknownMeanPrecision::MultivariateNormalUnknownMeanPrecision(const Vector& mu, const real tau, const real alpha, const Matrix& T) 
    : n_dim(mu.Size()),
      marginal_mean(n_dim),
      marginal(n_dim),
      mu_0(mu),
      tau_0(tau),
      mu_n(n_dim),
      alpha_0(alpha),
      T_0(T),
      bx_n(n_dim),
      M_2n(n_dim, n_dim)
{
    Reset();
}
void MultivariateNormalUnknownMeanPrecision::Reset()
{
    marginal_mean.setDegrees(alpha_0 - (real) n_dim + 1.0);
    marginal_mean.setLocation(mu_0);
    marginal_mean.setPrecision((alpha_0 - (real) n_dim + 1.0) * T_0.Inverse());
    n = 0;
    sum = 0.0;
    tau_n = tau_0;
    mu_n = mu_0;
    alpha_n = alpha_0;
    T_n = T_0;
    M_2n.Clear();
    bx_n.Clear();
}
MultivariateNormalUnknownMeanPrecision::~MultivariateNormalUnknownMeanPrecision()
{
}
Vector MultivariateNormalUnknownMeanPrecision::generate()
{
    Serror("Fix me!\n");
    return mu_n;
}
Vector MultivariateNormalUnknownMeanPrecision::generate() const
{
    Serror("Fix me!\n");
    return mu_n;
}

real MultivariateNormalUnknownMeanPrecision::Observe(const Vector& x)
{
    real p = marginal_mean_pdf(x);
    calculatePosterior(x);
    return p;
}


/// This is the probability of a particular Multi-Variate Nromal
real MultivariateNormalUnknownMeanPrecision::pdf(const Vector& mean) const
{
    Serror("Not implemented\n");
    return 0;
}



/** The marginal pdf of the observations.
    
    Instead of calculating the actual marginal:
    \f[
    \xi(x) = \int f(x \mid m, r) \, d\xi(m, r),
    \f]
    we calculate:
    \f[
    \xi(x) = f(x \mid E_\xi m, E_\xi r),
    \f]
    where \f$E_\xi m = \int m d\xi(m) \f$, \f$E_\xi r = \int r d\xi(r) \f$.
*/
real MultivariateNormalUnknownMeanPrecision::marginal_pdf(const Vector& x) const
{
    real p = marginal.pdf(x);
    if (std::isnan(p)) {
        marginal.Show();
    }
    return p;


}
real MultivariateNormalUnknownMeanPrecision::log_pdf(const Vector& x) const
{
    Serror("Not Implemented\n");
    return 0.0;
}


const Vector& MultivariateNormalUnknownMeanPrecision::getMean() const
{
    return mu_n;
}


/** Calculate the posterior parameter distribution given a new real observation.

    As opposed to the unimodal case, \f$x \in R^d\f$.
    We need to calculate the \f$d \times d\f$ matrix:
    \f[
    M_{2,n} = \sum_{i=1}^n (x_i - \bar{x}_n)^2
    \f]
    at each pass. However we use Knuth's algorithm to obtain:
    \f[
    M_{2,n} = M_{2,n-1} + (x_n - \bar{x}_{n})(x_n - \bar{x}_{n-1})'
    \f]

    
 */
void MultivariateNormalUnknownMeanPrecision::calculatePosterior(const Vector& x)
{
    n++;
    Vector bxn_prev = bx_n;
    real rn = (real) n;
    real irn = 1 / rn;
    bx_n = bx_n + (x - bx_n) * irn;
    for (int i=0; i<n_dim; ++i) {
        for (int j=0; j<n_dim; ++j) {
            M_2n(i,j) += (x[i] - bx_n[i])*(x[j] - bxn_prev[j]);
        }
    }
    Vector delta_mean = (bx_n - mu_0);
    
    real gamma = (rn - 1) * irn;
    mu_n = mu_n * gamma + x * (1 - gamma);
    tau_n += 1.0;
    alpha_n += 1.0;
    // T_n = T_0 + M_2n + delta_mean*delta_mean' * tau_0 * n / tau_n;
    Product(&delta_mean, &delta_mean, &T_n);
    T_n = T_0 + M_2n + T_n * tau_0 * n / tau_n;
    
    Matrix InvT = T_n.Inverse();
    marginal.setDegrees(alpha_n - (real) n_dim + 1.0);
    marginal.setLocation(mu_n);
    marginal.setPrecision(tau_n * (alpha_n - (real) n_dim + 1.0) * InvT);
    //    p_x_mr.setMean(mu_n);
    //    p_x_mr.setAccuracy(n * InvT);
}


void MultivariateNormalUnknownMeanPrecision::Show() 
{
    printf("# Norma-Wishart - ");
    printf("# Conditional normal mean: ");  mu_n.print(stdout);
    printf("# Conditional normal scaling: %f\n", tau_n);
    printf("# Wishart Degrees: %f\n", alpha_n);
    printf("# Wishart Precision: \n");  T_n.print(stdout);
    printf("# Normal-Wishart Predictive Marginal:\n");
    marginal.Show();
}
