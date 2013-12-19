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

#ifndef SPARSE_GAUSSIAN_PROCESS_H
#define SPARSE_GAUSSIAN_PROCESS_H

#include "Vector.h"
#include "Matrix.h"
#include "Distribution.h"
#include "NormalDistribution.h"

#include <vector>

/** Sparse Gaussian process. */
class SparseGaussianProcess
	{
	protected:
		Matrix X; ///< Samples
		Vector Y; ///< Output
		int N;  ///< Total number of basis funcitons (Dictionary size).

		Matrix Sigma_p;
		Matrix L;  ///< Cholesky Decomposition (L is an upper tringular matrix)
		Matrix inv_L;
		Matrix K;  ///< Kernel(Covariance) Matrix
		Matrix inv_K; ///< Inverse Covariance Matrix (Precision)
		/// Kernel hyperparameters.
		real noise_variance;	///< noise variance
		Vector scale_length;	///< lenght scale
		real sig_var;		    ///< signal variance 
		real v; ///<Dictionary threshold
		Vector scale_length_dic; ///< scale length for dictionary
		Vector alpha;
	public:
		SparseGaussianProcess(real noise_variance_,
							  Vector scale_length_,
							  real sig_var_,
							  real threshold_,
							  Vector scale_length_dic_);
		SparseGaussianProcess(Matrix& X_, 
							  Vector& Y_,
							  real noise_variance_,
							  real scale_length_,
							  real sig_var_,
							  real threshold_,
							  Vector scale_length_dic_);
		virtual ~SparseGaussianProcess();
		virtual int getNSamples() { return N; }
		virtual Vector generate();
		virtual real pdf(Vector& x, real y);
		virtual void Observe(const Matrix& x, const Vector& y);
		virtual void Observe(const std::vector<Vector>& x, const std::vector<real>& y);
		virtual void AddObservation(const Vector& x, const real& y);
		virtual void AddObservation(const std::vector<Vector>& x, const std::vector<real>& y);
		virtual void UpdateSparseGaussianProcess();
		virtual real GeneratePrediction(const Vector& x);
		virtual real GeneratePredictionKernel(const Vector& k);
		virtual void Prediction(Vector& x, real& mean, real& var);
		virtual real PredictiveMean(const Vector& x);
		virtual real PredictiveVariance(const Vector& x);
		virtual void Covariance();
		virtual Vector Kernel(const Vector& x);
		virtual Vector Kernel(const Vector& x, const Vector& scale);
		virtual void Clear();
	};

#endif


