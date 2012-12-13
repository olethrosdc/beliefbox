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

#include "BayesianMultivariateRegression.h"
#include "SpecialFunctions.h"

BayesianMultivariateRegression::BayesianMultivariateRegression(int m_, int d_, Matrix S0_, real N0_, real a_)
	:m(m_), d(d_), S0(S0_), N0(N0_), a(a_)
{
	N = 0;
	M = Matrix::Null(d,m);
	A = M;
	V = Matrix::Unity(d,d);
	K = Matrix::Unity(m,m)*a;
	Sxx = K;
	inv_Sxx = Sxx.Inverse();
	Syx = M*K;
	Syy = Syx * Transpose(M);
	Sy_x = Syy - (Syx*inv_Sxx)*Transpose(Syx);
}

void BayesianMultivariateRegression::AddElement(const Vector& y, const Vector& x)
{
	assert(y.Size() == d);
	assert(x.Size() == m);
	N++;
	Sxx = Sxx + OuterProduct(x,x);
	inv_Sxx = Sxx.Inverse();
	Syx = Syx + OuterProduct(y,x);
	Syy = Syy + OuterProduct(y,y);
	K = Sxx;
	M = Syx*inv_Sxx;
	Sy_x = Syy - M*Transpose(Syx);
}

/// Generate response matrix
Matrix BayesianMultivariateRegression::generate()
{
	iWishart iwishart(N + N0, Sy_x + S0, true);
	V = iwishart.generate();
	MultivariateNormal multivariate_normal(M.Vec(), Kron(V.Inverse(),K));
	Vector mean = multivariate_normal.generate();

	A.Vec(mean);

	return A;
}

void BayesianMultivariateRegression::generate(Matrix& MM, Matrix& VV)
{
	iWishart iwishart(N + N0, Sy_x + S0, true);
	VV = iwishart.generate();
	MultivariateNormal multivariate_normal(M.Vec(), Kron(VV.Inverse(),K));
	Vector mean = multivariate_normal.generate();
	MM.Vec(mean);
}

real BayesianMultivariateRegression::PredictiveDistribution(const Vector& x, const Vector& y)
{
	assert(x.Size() == m);
	
	real n = (N + N0 + 1.0) / 2.0;

    const Matrix mean = M;
	const Vector residual = y - mean*x;

	const Matrix inv_Sxx_new = (Sxx + OuterProduct(x,x)).Inverse();
	const Vector xSxx = inv_Sxx_new*x;
	real c = Product(x,xSxx);
	c = 1.0 - c;
	const Matrix r1 = (Sy_x.Inverse()*c);
	Vector r2 = r1*residual;
	
	return (Gamma(n) / Gamma(n - d/2.0)) * (1.0 / sqrt((Sy_x * (M_PI / c)).det()))*pow(Product(residual,r2) + 1.0 , -n);
}

void BayesianMultivariateRegression::Reset()
{
	N = 0;
	M = Matrix::Null(d,m);
	A = M;
	V = Matrix::Unity(d,d);
	K = Matrix::Unity(m,m)*a;
	Sxx = K;
	inv_Sxx = Sxx.Inverse();
	Syx = M*K;
	Syy = Syx * Transpose(M);
	Sy_x = Syy - (Syx*inv_Sxx)*Transpose(Syx);
}
