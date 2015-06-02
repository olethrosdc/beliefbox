// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include "GaussianProcessBandit.h"
 
GaussianProcessBandit::GaussianProcessBandit(Vector gp_scale_length_, 
	Vector gp_scale_length_dic_, real gp_noise_var_, real gp_sig_var_, 
	real gp_threshold_, real delta_) 
	: gp(gp_noise_var_, gp_scale_length_, gp_sig_var_, 
		 gp_threshold_, gp_scale_length_dic_)
{	
	reset();
	delta = delta_;
}
 
GaussianProcessBandit::~GaussianProcessBandit() 
{ 
}
real GaussianProcessBandit::beta() 
{ 
	return 2 * std::log(d * trial * trial * M_PI * M_PI / (6 * delta));
}

int GaussianProcessBandit::predict(const Matrix& X_, Vector& predictionVector) 
{
	int N = X_.Rows();
	predictionVector.Resize(N);
	
	real value;
	real mean;
	real var;
	Vector row;
	
	for(int i=0;i<N;i++) {
		row = X_.getRow(i);
		gp.Prediction(row,mean,var);
		value = mean + sqrt(beta()) * var; //UCB
		predictionVector(i) = value;
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
void GaussianProcessBandit::addObservation(const Matrix& X_, const Vector& Y_)
{
	int N = X_.Rows();
	d = X_.Columns();
	std::vector<Vector> X(N);
	std::vector<real> Y(N);
	for(int i=0;i<N;i++) {
		X[i] = X_.getRow(i);
		Y[i] = Y_(i);
	}
	gp.AddObservation(X,Y);
}
void GaussianProcessBandit::addObservation(const Vector& x_, const real& y_)
{
	d = x_.Size();
	std::vector<Vector> X(1);
	std::vector<real> Y(1);
	X[0] = x_;
	Y[0] = y_;
	gp.AddObservation(X,Y);
}
void GaussianProcessBandit::reset() 
{
	gp.Clear();
	trial = 1;
}