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
#include "Student.h"

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

Vector BayesianMultivariateRegression::generate(const Vector& x)
{	
	Matrix mean = generate();

	return mean*x;
}

void BayesianMultivariateRegression::generate(Matrix& MM, Matrix& VV)
{
	iWishart iwishart(N + N0, Sy_x + S0, true);
	VV = iwishart.generate();

	MultivariateNormal multivariate_normal(M.Vec(), Kron(VV.Inverse(),K));

	Vector mean = multivariate_normal.generate();
	MM.Vec(mean);
}

real BayesianMultivariateRegression::Posterior(const Vector& x, const Vector& y)
{
	iWishart iwishart(N + N0, Sy_x + S0, true);
	Matrix VV = iwishart.generate();
	Matrix MM = Matrix::Null(d,m);
	
	MultivariateNormal multivariate_normal(M.Vec(), Kron(VV.Inverse(),K));
	Vector mean = multivariate_normal.generate();
	MM.Vec(mean);
	
	MultivariateNormal m_normal((MM*x),VV.Inverse());

	return m_normal.pdf(y);
	//real result = LogPredictiveDistribution(x,y);
//	printf("Log Prediction = %f and Prediction = %f\n",result,exp(result));
//	return exp(result);
}

real BayesianMultivariateRegression::PredictiveDistribution(const Vector& x, const Vector& y)
{
	assert(x.Size() == m);
	real n = (N + N0 + 1.0) / 2.0;

    const Matrix mean = M;
	const Vector residual = y - mean*x;

	Matrix xxx = OuterProduct(x,x);
	const Matrix inv_Sxx_new = (Sxx + xxx).Inverse();
	const Vector xSxx = inv_Sxx_new*x;
	real c = Product(x,xSxx);
	c = 1.0 - c;

	Matrix r1 = Sy_x.Inverse();
	r1 = r1*c;
	Vector r2 = r1*residual;
	real r3 = Product(residual, r2);
	real r4 = (Sy_x * (M_PI / c)).det();
//	real r5 = Gamma(n);
	
//	return (Gamma(n) / Gamma(n - d/2.0)) * (pow(r4,-0.5))*pow((r3+ 1.0) , -n);
	return (pow(r4,-0.5))*pow((r3+ 1.0) , -n);	
}

real BayesianMultivariateRegression::LogPredictiveDistribution(const Vector& x, const Vector& y) {
	//Degree 
	real degrees = (real)(N + N0 + 1);
	//Mean
	Vector mu = M*x;
	printf("Output\n");
	y.print(stdout);
	printf("Mean\n");
	mu.print(stdout);
	//Precision 
	Matrix t = Sxx + OuterProduct(x,x);
	t = t.Inverse_LU();
	Vector Temp = t * x;	
	real c	= Product(x,Temp);
	c		= 1.0 - c;
	
//	c =1.0;
	printf("Precision\n");
	Matrix precision1 = (Sy_x + S0).Inverse();
	precision1.print(stdout);
//	printf("C = %f\n",c);
	Matrix precision = (Sy_x + S0).Inverse() * (c);
	Matrix variance  = ((Sy_x + S0)*((1.0/c)));
	
	printf("Precision\n");
	precision.print(stdout);
	printf("Degree => %f\n",degrees);
//	Student s(degrees, mu, precision);
//	return s.pdf(y);
	
	Vector delta = y - mu;
	real g = 1 + Mahalanobis2(delta, precision, delta)/degrees;
//	real g = 1 + Mahalanobis2(delta, precision, delta);

	real d = (real)y.Size();
	
	real log_c =  logGamma(0.5 *(degrees + d) )
				- 0.5 * log(variance.det())
				- logGamma(0.5 * degrees)
				- (0.5 * d) * log(degrees*M_PI);
	real log_p = log_c - (0.5*degrees)*log(g);
	//real log_c = logGamma(0.5 * degrees) - 0.5 * log(variance.det()) - logGamma(0.5* degrees - 0.5*d);
//	real log_p = log_c - (0.5*degrees)*g;
	return log_p;
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
