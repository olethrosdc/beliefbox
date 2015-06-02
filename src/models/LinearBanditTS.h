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
  * Note: This is a non-standard way of doing Thompson sampling
  *
  * Based on the work in
  * Thompson Sampling for Contextual Bandits with Linear Payoffs
  * by Agrawal and Goyal
  */
 
#ifndef LINEARTS_CONTEXT_BANDIT_H
#define LINEARTS_CONTEXT_BANDIT_H

#include "real.h"
#include "Matrix.h"
#include "Vector.h"

class LinearBanditTS
{
private:
	Matrix B;
	Vector f;
	Vector meanHat;
	
	real delta;
	real epsilon;
	real R;
	
	bool empty;
	int trial;
	int d;
	real V();
	void initializeBandit();
	void updateBandit(int index, const Vector& x_, const real& y_);
	
public:
	LinearBanditTS(real R_, real delta_, real epsilon_);
	~LinearBanditTS();
	int predict(const Matrix& X_, Vector& predictionVector); 
		//predicts the next index to test;
	void addObservation(const Matrix& X_, const Vector& Y_); 
	void addObservation(const Vector& x_, const real& y_);
	void reset();
};

#endif