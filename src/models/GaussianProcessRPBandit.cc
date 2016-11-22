// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include "GaussianProcessRPBandit.h"
#include "Random.h"
 
GaussianProcessRPBandit::GaussianProcessRPBandit(Vector gp_scale_length_, 
	Vector gp_scale_length_dic_, real gp_noise_var_, real gp_sig_var_, 
	real gp_threshold_, real delta_) 
	: gp(gp_noise_var_, gp_scale_length_, gp_sig_var_, 
		 gp_threshold_, gp_scale_length_dic_)
{	
	reset();
	delta = delta_;
	gp_scale_length = gp_scale_length_;
	gp_scale_length_dic = gp_scale_length_dic_;
	gp_noise_var = gp_noise_var_;
	gp_sig_var = gp_sig_var_;
	gp_threshold = gp_threshold_;

	X = Matrix();
	Y = Vector();
}

GaussianProcessRPBandit::~GaussianProcessRPBandit() 
{ 
}
real GaussianProcessRPBandit::beta(int depth) 
{ 
	int T = trial + depth;
	return 2 * std::log(d * T * T * M_PI * M_PI / (6 * delta));
}

int GaussianProcessRPBandit::predict(const Matrix& X_, Vector& predictionVector, int projections, int depth) 
{
	int N = X_.Rows();
	int n = X.Rows();
	predictionVector = predictionVector.Null(N);

	SparseGaussianProcess tmpGP(gp_noise_var,gp_scale_length,gp_sig_var,gp_threshold,gp_scale_length_dic);
	
	real value;
	real mean;
	real var;
	Vector row;
	
	for(int i=0;i<projections;i++) {
		tmpGP.Observe(X,Y);
		for(int j=0;j<depth;j++) {
			int index = urandom(0,X_.Rows());
			tmpGP.Prediction(X_.getRow(index),mean,var);
			value = mean + sqrt(beta(j)) * var; //UCB

			std::vector<Vector> x_(1);
			std::vector<real> y_(1);
			x_[0] = X_.getRow(index);
			y_[0] = value;
			tmpGP.AddObservation(x_,y_);
		}
		for(int j=0;j<N;j++) {
			tmpGP.Prediction(X_.getRow(j),mean,var);
			value = mean + sqrt(beta(depth)) * var;
			if(value > predictionVector(j))
				predictionVector(j) = value;
		}
	}
	int bestIndex = 0;
	real maxVal = predictionVector(0);
	for(int i=0;i<N;i++) {
		if(predictionVector(i) > maxVal)
		{
			maxVal = predictionVector(i);
			bestIndex = i;
		}
	}
	trial++;
	return bestIndex;
}
void GaussianProcessRPBandit::addObservation(const Matrix& X_, const Vector& Y_)
{
	int N = X_.Rows();
	d = X_.Columns();
	if(empty) 
	{
		X.Resize(N,d);
	}
	std::vector<Vector> X__(N);
	std::vector<real> Y__(N);
	for(int i=0;i<N;i++) {
		X__[i] = X_.getRow(i);
		Y__[i] = Y_(i);

		X = X.AddRow(X_.getRow(i));
		Y.AddElement(Y_(i));
	}
	gp.AddObservation(X__,Y__);
	empty = false;
}
void GaussianProcessRPBandit::addObservation(const Vector& x_, const real& y_)
{
	if(empty) {
		X.Resize(1,x_.Size());
		X.setRow(0,x_);
	}
	else {
		X = X.AddRow(x_);
	}
	d = x_.Size();
	std::vector<Vector> X__(1);
	std::vector<real> Y__(1);
	X__[0] = x_;
	Y__[0] = y_;
	Y.AddElement(y_);
	gp.AddObservation(X__,Y__);
	empty = false;
}
void GaussianProcessRPBandit::reset() 
{
	gp.Clear();
	X = Matrix();
	Y = Vector();
	trial = 1;
	empty = true;
}