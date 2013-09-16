/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GAUSSIAN_PROCESS_H
#define GAUSSIAN_PROCESS_H

#include "Vector.h"
#include "Matrix.h"
#include "Distribution.h"
#include "NormalDistribution.h"

#include <vector>

/** Gaussian process. 
    
    This is a {\em conditional} distribution.
 */
class GaussianProcess
{
protected:
	Matrix X; ///< Samples
	Vector Y; ///< Output
	int N;  ///< Total number of samples
	Vector alpha; 
    Matrix Sigma_p;
    Matrix Accuracy;
    Matrix A;
	Matrix L;  ///< Cholesky Decomposition (L is an upper tringular matrix).
	Matrix inv_L;
	Matrix K;  ///< Kernel(Covariance) Matrix.
	/// Kernel hyperparameters.
    real noise_variance;	///< noise variance
	Vector scale_length;	///< lenght scale
	real sig_var;		    ///< signal variance 
    Matrix X2; ///< observation co-variance
    Vector mean;
    Matrix covariance;
public:
    GaussianProcess(Matrix& Sigma_p_,
                    real noise_variance_);
	GaussianProcess(real noise_variance_,
					Vector scale_length_,
					real hyp_u_);
	GaussianProcess(Matrix& X_, 
					Vector& Y_,
					real noise_variance_,
					real scale_length_,
					real hyp_u_);
    virtual ~GaussianProcess();
    virtual Vector generate();
    virtual real pdf(Vector& x, real y);
    virtual void Observe(Vector& x, real y);
    virtual void Observe(Matrix& x, Vector& y);
	virtual void UpdateGaussianProcess();
	virtual real GeneratePrediction(const Vector& x);
	virtual void Prediction(Vector& x, real& mean, real& var);
	virtual real PredictiveMean(const Vector& x);
	virtual real PredictiveVariance(const Vector& x);
	virtual void Covariance();
	virtual Matrix CovarianceDerivatives(int p);
	virtual Vector Kernel(const Vector& x);
	virtual real LogLikelihood();
};

#endif

