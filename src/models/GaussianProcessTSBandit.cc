// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include "GaussianProcessTSBandit.h"
#include "Random.h"
 
GaussianProcessTSBandit::GaussianProcessTSBandit(Vector gp_scale_length_, 
	Vector gp_scale_length_dic_, real gp_noise_var_, real gp_sig_var_, 
	real gp_threshold_, real delta_) 
	: gp(gp_noise_var_, gp_scale_length_, gp_sig_var_, 
		 gp_threshold_, gp_scale_length_dic_),
	  tmpGP(gp_noise_var_, gp_scale_length_, gp_sig_var_, 
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

GaussianProcessTSBandit::~GaussianProcessTSBandit() 
{ 
}
real GaussianProcessTSBandit::beta(int depth) 
{ 
	return 2 * std::log(d * depth * depth * M_PI * M_PI / (6 * delta));
}
real GaussianProcessTSBandit::descend(const Matrix& X_, const Matrix& learnedMatrix, const Vector& learnedVector, int index, int depth, int horizon, real probability) {
	real mean;
	real var;
	real value;
	tmpGP.Observe(learnedMatrix, learnedVector);
	NormalDistribution normal;
	if(depth==horizon) {
		real maxVal = 0;
		for(int i=0;i<X_.Rows();i++) {
			tmpGP.Prediction(X_.getRow(i),mean,var);
			normal = NormalDistribution(mean,sqrt(var));
			value = normal.generate();
			//value = mean + sqrt(beta(X.Rows())) * var; //UCB
			if(value > maxVal) 
				maxVal = value;
		}
		return maxVal;
	}
	else {
		int N = X_.Rows();
		int count = 0;
		real QVal = 0;
		for(int i=0;i<N;i++) {
			if(probability > urandom()) { //descend
				count++;
				tmpGP.Prediction(X_.getRow(i),mean,var);
				normal = NormalDistribution(mean,sqrt(var));
				value = normal.generate();
				//value = mean + sqrt(beta(X.Rows())) * var; //UCB
				Matrix tmpMatrix = Matrix(learnedMatrix);
				Vector tmpVector = Vector(learnedVector);
				tmpMatrix = tmpMatrix.AddRow(X_.getRow(i));
				tmpVector.AddElement(value);
				real tmp = descend(X_,tmpMatrix,tmpVector,i,depth+1,horizon,probability);
				if(tmp>QVal)
					QVal = tmp;
			}
		}
		if(count==0) {
			for(int i=0;i<X_.Rows();i++) {
				tmpGP.Prediction(X_.getRow(i),mean,var);
				normal = NormalDistribution(mean,sqrt(var));
				value = normal.generate();
				//value = mean + sqrt(beta(X.Rows())) * var; //UCB
				if(value > QVal) 
					QVal = value;
			}
		}
		return QVal;
	}
}
int GaussianProcessTSBandit::predict(const Matrix& X_, Vector& predictionVector, int projections, int horizon, real probability) 
{
	int N = X_.Rows();
	int n = X.Rows();
	predictionVector = predictionVector.Null(N);

	SparseGaussianProcess tmpGP(gp_noise_var,gp_scale_length,gp_sig_var,gp_threshold,gp_scale_length_dic);
	
	real value;
	real mean;
	real var;
	Vector row;
	
	for(int i=0;i<N;i++) {
		tmpGP.Observe(X,Y);
		if(probability > urandom()) { //descend
			predictionVector(i) = descend(X_,X,Y,i,0,horizon,probability);
		}
		else {		
			tmpGP.Prediction(X_.getRow(i),mean,var);
			NormalDistribution normal = NormalDistribution(mean,sqrt(var));
			value = normal.generate();
			predictionVector(i) = value;
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
	return bestIndex;
}
void GaussianProcessTSBandit::addObservation(const Matrix& X_, const Vector& Y_)
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
void GaussianProcessTSBandit::addObservation(const Vector& x_, const real& y_)
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
void GaussianProcessTSBandit::reset() 
{
	gp.Clear();
	X = Matrix();
	Y = Vector();
	empty = true;
}