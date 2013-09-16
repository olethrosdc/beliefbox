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

GaussianProcess::GaussianProcess(real noise_variance_,
								 Vector scale_length_,
								 real sig_var_)
	: noise_variance(noise_variance_),
	  scale_length(scale_length_),
	  sig_var(sig_var_)
{
}

GaussianProcess::GaussianProcess(Matrix& X_, 
								 Vector& Y_,
								 real noise_variance_,
								 real scale_length_,
								 real sig_var_)
	: X(X_),
	  Y(Y_),
	  noise_variance(noise_variance_),
	  scale_length(scale_length_),
	  sig_var(sig_var_)
{
	N = X.Rows();
	UpdateGaussianProcess();
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
    MatrixProduct(&x, &x, &V);
    X2 += V;

    A = X2 / noise_variance + Accuracy;
    Matrix inv_A = A.Inverse();
    //mean = inv_A * X;
}


/// This implements functon space view of a GP
///
/// X has a number of columns equal to the amount of data.
//void GaussianProcess::Observe(Matrix& X, Vector& y)
//{
//    int N = X.Columns();
//    Matrix K(N, N);
//    for (int i=0; i<N; ++i) {
//        Vector S = K.getColumn(i);
//        for (int j=0; j<N; ++j) {
//            real delta = (S - K.getColumn(j)).SquareNorm();
//            K(i,j) = exp(-0.5 * delta);
//            if (i==j) {
//                K(i,j) += noise_variance;
//            }
//        }
//    }
//	
//    Matrix L = K.Cholesky();
//}
void GaussianProcess::Observe(Matrix& X_, Vector& Y_)
{
	X = X_;
	Y = Y_;
    N = X.Rows();
	UpdateGaussianProcess();
}

void GaussianProcess::UpdateGaussianProcess()
{
	Covariance();
	L = K.Cholesky();
	inv_L = L.Inverse();
	alpha = inv_L*(Transpose(inv_L)*Y); 
}

real GaussianProcess::GeneratePrediction(const Vector& x)
{
	Vector k  = Kernel(x);
	real mean = PredictiveMean(k);
	real var  = PredictiveVariance(k);
	return mean;
	//NormalDistribution N(mean, sqrt(var));
//	return N.generate();
}

void GaussianProcess::Prediction(Vector& x, real& mean, real& var)
{
	Vector k = Kernel(x);
	mean     = PredictiveMean(k);
	var	     = PredictiveVariance(k);
}

real GaussianProcess::PredictiveMean(const Vector& k)
{
	real mean = Product(k, alpha); ///(mean = k'*alpha) 
	return mean;
}

real GaussianProcess::PredictiveVariance(const Vector& k)
{
	Vector v = Transpose(inv_L)*k;
	real var = sig_var*sig_var - Product(v,v);
	return var;
}

/// Covariance function estimation
void GaussianProcess::Covariance()
{
	/// The covariance matrix is symmetric
	K = Matrix::Matrix(N,N);
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
//	K.print(stdout);
}

/// The following function calculates the partial derivatives of the Covariance function w.r.t hyperparameters
/// If input argument equal to 1, the derivative of the scale_lenght 'hyperparameter' is calculated.
/// Otherwise, the derivative of the signal variance 'hyperparameter' is calculated.
Matrix GaussianProcess::CovarianceDerivatives(int p)
{
	Matrix KK(N,N);
	Matrix KE(N,N);
	for(int i=0; i<N; ++i) { 
		Vector S = X.getRow(i);
		for(int j = i; j < N; ++j) {
			real delta = ((S - X.getRow(j))/scale_length).SquareNorm();
			KK(i,j) = delta;
			KK(j,i) = KK(i,j);
			KE(i,j) = exp(-0.5*delta);
			KE(j,i) = KE(i,j);
		}
	}
	
	Matrix DK(N,N);
	if(p == 1)  
	{
		DK = (sig_var*sig_var)*KE.Multiple(KK);
	}
	else if(p == 2) {
		DK = 2*sig_var*KE;
	}
	return DK;
}

/// Covariance funtion de
Vector GaussianProcess::Kernel(const Vector& x)
{
	//int N = X.Rows();
	Vector k(N);
	for(int i=0; i<N; ++i) { 
		real delta = ((x - X.getRow(i))/scale_length).SquareNorm();
		delta = sig_var*sig_var*exp(-0.5*delta);
		k(i) = delta; 
	}
	return k;
}

/// Log marginal likelihood computation
real GaussianProcess::LogLikelihood()
{
	real slogL = 0.0;
	for(int i=0; i<N; ++i) {
		slogL += log(L(i,i));
	}
	real LogLik = -0.5*Product(Y, alpha) - slogL - (0.5*N)*log(2*M_PI);
	return LogLik;
}

