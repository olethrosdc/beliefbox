/* -*- Mode: C++; -*- */
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SparseGaussianProcess.h"

/// Create a new GP with observations
SparseGaussianProcess::SparseGaussianProcess(real noise_variance_,
											 Vector scale_length_,
											 real sig_var_,
											 real threshold_,
											 Vector scale_length_dic_)
	: noise_variance(noise_variance_),
	scale_length(scale_length_),
	sig_var(sig_var_), 
	v(threshold_),
	scale_length_dic(scale_length_dic_)
{
	N = 0;
}

SparseGaussianProcess::SparseGaussianProcess(Matrix& X_, 
											 Vector& Y_,
											 real noise_variance_,
											 real scale_length_,
											 real sig_var_,
											 real threshold_,
											 Vector scale_length_dic_)
	: X(X_),
	Y(Y_),
	noise_variance(noise_variance_),
	scale_length(scale_length_),
	sig_var(sig_var_),
	v(threshold_),
	scale_length_dic(scale_length_dic_)
{
	Observe(X_, Y_);
}

SparseGaussianProcess::~SparseGaussianProcess()
{
}

Vector SparseGaussianProcess::generate()
{
    Serror ("Not implemented\n");
    exit(-1);
    //return Vector();
}    

real SparseGaussianProcess::pdf(const Vector& x, real y)
{
    return 0.0;
}

void SparseGaussianProcess::Observe(const Matrix& x, const Vector& y)
{
	N = 0;
	for(int i = 0; i<x.Rows(); i++) {
		AddObservation(x.getRow(i),y(i));
	}
	printf("The dictionary contains N = %d basis functions\n",N);
	UpdateSparseGaussianProcess();
}

void SparseGaussianProcess::Observe(const std::vector<Vector>& x, const std::vector<real>& y)
{
	N = 0;
	for(int i = 0; i<(int)x.size(); ++i) {
		AddObservation(x[i], y[i]);
	}
	printf("The dictionary contains N = %d basis functions\n",N);
	UpdateSparseGaussianProcess();
}

void SparseGaussianProcess::AddObservation(const Vector& x, const real& y)
{
	Vector k;
	Vector iKk;
	real mu;
	if(N == 0) {
		inv_K = Matrix::Unity(1,1);
		inv_K(0,0) = 1.0/(sig_var*sig_var);
		N = 1;
		X = Matrix(1,x.Size());
		Y = Vector(1);
		X.setRow(0,x);
		Y(0) = y;
	} else {
		k = Kernel(x, scale_length_dic);
		iKk = inv_K*k;
		real rr =  Product(iKk, k);
		mu = (1.0 - rr);
		if(mu > v) {
			inv_K = (mu*inv_K) + OuterProduct(iKk,iKk);
			iKk = iKk*(-1.0);
			inv_K = inv_K.AddRow(iKk);
			iKk.AddElement(1.0);
			inv_K = inv_K.AddColumn(iKk);
			Matrix kk = inv_K*(1/mu);
			inv_K= kk;
			N++;
			X = X.AddRow(x);
			Y.AddElement(y);
		}
	}
}

void SparseGaussianProcess::AddObservation(const std::vector<Vector>& x, const std::vector<real>& y)
{
	for(int i = 0; i<(int)x.size(); ++i) {
		AddObservation(x[i], y[i]);
	}
	printf("The dictionary contains N = %d basis functions\n",N);
	UpdateSparseGaussianProcess();
}

void SparseGaussianProcess::UpdateSparseGaussianProcess()
{
	Covariance();
	L = K.Cholesky();
	inv_L = L.Inverse();
//	alpha = inv_L*(Transpose(inv_L)*Y); 
	alpha = K.Inverse()*Y; 
}

real SparseGaussianProcess::GeneratePrediction(const Vector& x)
{
	if(N > 0) {
		Vector k  = Kernel(x);
		real mean = PredictiveMean(k);
		real var  = PredictiveVariance(k);
//		return mean;
		NormalDistribution N(mean, sqrt(var));
		return N.generate();
	}else {
		NormalDistribution N(0, 1);
		return N.generate();
	}
}

real SparseGaussianProcess::GeneratePredictionKernel(const Vector& k)
{
	if(N > 0) {
		real mean = PredictiveMean(k);
		real var  = PredictiveVariance(k);

		NormalDistribution N(mean, sqrt(var));
		return N.generate();
	}else {
		NormalDistribution N(0, 1);
		return N.generate();
	}
}

void SparseGaussianProcess::Prediction(const Vector& x, real& mean, real& var)
{
	Vector k = Kernel(x);
	mean     = PredictiveMean(k);
	var	     = PredictiveVariance(k);
}

real SparseGaussianProcess::PredictiveMean(const Vector& k)
{
	real mean = Product(k, alpha); ///(mean = k'*alpha) 
	return mean;
}

real SparseGaussianProcess::PredictiveVariance(const Vector& k)
{
	Vector iLk = Transpose(inv_L)*k;
	real var = sig_var*sig_var - Product(iLk,iLk);
	return var;
}

/// Covariance function estimation
void SparseGaussianProcess::Covariance()
{
	/// The covariance matrix is symmetric
	K = Matrix(N,N);
	for(int i=0; i<N; ++i) { 
		Vector S = X.getRow(i);
		for(int j = i; j < N; ++j) {
			real delta = ((S - X.getRow(j))/scale_length).SquareNorm();
			delta = sig_var*sig_var*exp(-0.5*delta);
			if(i == j) {
				K(i,j) = delta + noise_variance*noise_variance;
			}else {
				K(i,j) = delta;
				K(j,i) = delta;
			}
		}
	}
}

/// Kernel function
Vector SparseGaussianProcess::Kernel(const Vector& x)
{
	Vector k(N);
	real sig_noise = sig_var*sig_var;
	for(int i=0; i<N; ++i) { 
		real delta = ((x - X.getRow(i))/scale_length).SquareNorm();
		delta = sig_noise*exp(-0.5*delta);
		k(i) = delta; 
	}
	return k;
}

/// Kernel function
Vector SparseGaussianProcess::Kernel(const Vector& x, const Vector& scale)
{
	Vector k(N);
	real sig_noise = sig_var*sig_var;
	for(int i=0; i<N; ++i) { 
		real delta = ((x - X.getRow(i))/scale).SquareNorm();
		k(i) = sig_noise*exp(-0.5*delta);
	}
	return k;
}

void SparseGaussianProcess::Clear()
{
	N = 0;
}
