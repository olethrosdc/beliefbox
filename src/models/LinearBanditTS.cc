// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
 #include "LinearBanditTS.h"
 #include "MultivariateNormal.h"
 
 LinearBanditTS::LinearBanditTS(real R_, real delta_, real epsilon_)
 {
	R = R_;
	delta = delta_;
	epsilon = epsilon_;
	
	reset();
 }
 LinearBanditTS::~LinearBanditTS() 
 {
 }
 int LinearBanditTS::predict(const Matrix& X_, Vector& predictionVector) 
 {
	int N = X_.Rows();
	predictionVector.Resize(N);
	 
	MultivariateNormal mvn = MultivariateNormal(meanHat,
												pow(V(),2)*B.Inverse());
	Vector meanSample = mvn.generate();
	real v;
	for(int i=0;i<N;i++) {
		predictionVector(i) = Product(X_.getRow(i),meanSample);
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
 void LinearBanditTS::addObservation(const Matrix& X_, const Vector& Y_)
 {
	int N = X_.Rows();
	d = X_.Columns();
	
	if(empty) 
		initializeBandit();
	else {
		for(int i=0;i<N;i++) {
			updateBandit(i, X_.getRow(i), Y_(i));
		}
	}
	empty=false;
 }
 void LinearBanditTS::addObservation(const Vector& x_, const real& y_)
 {
	d = x_.Size();
	
	if(empty) 
		initializeBandit();
	else {
		updateBandit(0, x_, y_);
	}
	empty=false;
 }
 real LinearBanditTS::V() {
	 return R * sqrt(24 / epsilon * d * log(1 / delta));
 }
 void LinearBanditTS::initializeBandit()
 {
	meanHat = meanHat.Null(d);
	f = f.Null(d);
	B = B.Null(d,d);
	for(int i=0;i<d;i++) {
		B(i,i) = 1;
	}
 }
 void LinearBanditTS::updateBandit(int index, const Vector& x_, const real& y_) 
 {
	Matrix tmp = tmp.Null(d,d);
	MatrixProduct(x_,x_,tmp);
	B += tmp;
	f += x_ * y_;
	meanHat = B.Inverse() * f;
 }
 void LinearBanditTS::reset()
 {
	B = Matrix();
	f = Vector();
	meanHat = Vector();
	trial = 1;
	empty = true;
 }