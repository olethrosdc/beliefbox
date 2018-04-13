/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
// The specific class is based on the Thomas P. Minka paper "Bayesian Linear
// Regression", Technical Report, February 20, 2001.

#include "BayesianMultivariateRegression.h"
#include "SpecialFunctions.h"
#include "Student.h"
#include "gsl/gsl_sf_psi.h"

BayesianMultivariateRegression::BayesianMultivariateRegression(int m_, int d_, Matrix S0_, real N0_, real a_, bool ThompsonSampling_)
:m(m_), d(d_), S0(S0_), N0(N0_), a(a_), ThompsonSampling(ThompsonSampling_)
{
	N = 0.0;
	M = Matrix::Null(d,m);
	A = M;
	V = Matrix::Unity(d,d);
	K = Matrix::Unity(m,m)*a;
	Sxx = K;
	Syx = M*K;	
	Syy = Syx * Transpose(M);	
}

// y ~ N(Sx, V)
void BayesianMultivariateRegression::AddElement(const Vector& y, const Vector& x)
{
	assert(y.Size() == d);
	assert(x.Size() == m);
	N = N + 1.0;
	Sxx = Sxx + OuterProduct(x,x); // Sxx = X*X'
	Syx = Syx + OuterProduct(y,x); // Syx = Y*X' (Eq. 21)
	Syy = Syy + OuterProduct(y,y); // Syy = Y*Y' (Eq. 22)
	inv_Sxx = Sxx.Inverse_LU();  
	M = Syx*inv_Sxx;	
	Sy_x = Syy - M*Transpose(Syx); //Sy|x = Syy - Syx*Sxx^{-1}*Syx'	(Eq. 23)
}

/// Generate response matrix
Matrix BayesianMultivariateRegression::generate()
{
	Matrix S = A;
	if(N > 0) {
		iWishart iwishart(N + N0, Sy_x + S0, true); // Eq. 51
		V = iwishart.generate();
		MultivariateNormal multivariate_normal(M.Vec(), Kron(Sxx, V.Inverse())); // Eq. 10
		Vector mean = multivariate_normal.generate();
		S.Vec(mean);
	}	

	return S;
}

Vector BayesianMultivariateRegression::generate(const Vector& x) const
{	
	Vector new_state = A*x;
	return new_state;
}

void BayesianMultivariateRegression::generate(Matrix& MM, Matrix& VV)
{
	iWishart iwishart(N + N0, Sy_x + S0, true); // Eq. 51
	VV = iwishart.generate();
	
	MultivariateNormal multivariate_normal(M.Vec(), Kron(inv_Sxx,VV).Inverse_LU()); // Eq. 10
	
	Vector mean = multivariate_normal.generate();
	MM.Vec(mean);
}

/** Return the marginal posterior probability of y */
real BayesianMultivariateRegression::Posterior(const Vector& x, const Vector& y)
{
	Vector xx = (inv_Sxx*x);
	real c = 1.0 + Product(x,xx);
	V = ((Sy_x + S0)*c);
	Vector mean = M*x;
	Student st((N + N0 + 1.0), mean, V.Inverse_LU());
	
	return st.pdf(y);
}

void BayesianMultivariateRegression::Select()
{
	if(ThompsonSampling) {
		A = generate();
	} else {
		A = M;
	}
}

void BayesianMultivariateRegression::Reset()
{
	N = 0;
	M = Matrix::Null(d,m);
	A = M;
	V = Matrix::Unity(d,d);
	K = Matrix::Unity(m,m)*a;
	Sxx = K;
	Syx = M*K;	
	Syy = Syx * Transpose(M);	
}
