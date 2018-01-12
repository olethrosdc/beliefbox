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

#ifndef BAYESIAN_M_REGRESSION_H
#define BAYESIAN_M_REGRESSION_H

#include "BasisSet.h"
#include "Matrix.h"
#include "iWishart.h"
#include "MultivariateNormal.h"

///** Bayesian Multivariate Linear Regression **/
class BayesianMultivariateRegression
	{
	protected:
		int m; ///< Length of input vector
		int d; ///< Length of output vector
		real N; ///< Total number of samples
		Matrix A;
		Matrix M; ///< Mean
		Matrix V; ///< Rows covariance matrix
		Matrix K; ///< Columns covariance matrix
		Matrix S0; 
		real N0;
		real a;
		bool ThompsonSampling;
		///< Sufficient statistics
		Matrix Sxx;  ///< Matrix S_{xx}
		Matrix inv_Sxx;
		Matrix Syy;  ///< Matrix S_{yy}
		Matrix Syx;  ///< Matrix S_{yx}
		Matrix Sy_x; ///< Matrix S_{y|x} 
	public:
		BayesianMultivariateRegression(int m_ = 1, int d_ = 1, Matrix S0_ = Matrix::Unity(1,1), real N0_ = 1.0, real a_ = 1.0, bool ThompsonSampling_ = true);
		void AddElement(const Vector& y, const Vector& x);
		Matrix generate();
		Vector generate(const Vector& x) const;
		void generate(Matrix& Mean, Matrix& Covariance);
		Matrix getSxx() { return Sxx; }
		real Posterior(const Vector& x, const Vector& y);
		void Sampling(const bool ThompsonSampling_) { ThompsonSampling = ThompsonSampling_;}
		void Select();
		void Reset();
	};

#endif

