// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
 /*
  * Based on the work in
  * Gaussian Process Optimization in the Bandit Setting: 
		No Regret and Experimental Design
  * by Srinivas et al.
  */
 
#ifndef GAUSSIAN_PROCESS_TS_BANDIT_H
#define GAUSSIAN_PROCESS_TS_BANDIT_H

#include "SparseGaussianProcess.h"
#include "real.h"
#include "Matrix.h"
#include "Vector.h"

class GaussianProcessTSBandit
{
private: 
	SparseGaussianProcess gp;
	SparseGaussianProcess tmpGP;
	Vector gp_scale_length; //scaling for kernel
	Vector gp_scale_length_dic; //scaling for dictionary
	real gp_noise_var;
	real gp_sig_var;
	real gp_threshold; 

	Matrix X;
	Vector Y;
	bool empty;
	real delta;
	
	int d;
	real beta(int depth);
	real descend(const Matrix& X_, const Matrix& learnedMatrix, const Vector& learnedVector, int index, int depth, int horizon, real probability);
	
public:
	GaussianProcessTSBandit(Vector gp_scale_length, Vector gp_scale_length_dic, 
		real gp_noise_var, real gp_sig_var, real gp_threshold, real delta);
	~GaussianProcessTSBandit();
	int predict(const Matrix& X_, Vector& predictionVector, int projections, int horizon, real probability); 
		//predicts the next index to test;
	void addObservation(const Matrix& X_, const Vector& Y_); //learns GP
	void addObservation(const Vector& x_, const real& y_);
	void reset(); //generates a fresh GP
};

#endif
