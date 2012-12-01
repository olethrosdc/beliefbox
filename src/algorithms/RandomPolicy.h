/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_POLICY_H
#define RANDOM_POLICY_H

#include "AbstractPolicy.h"
#include "RandomNumberGenerator.h"

class RandomPolicy  : public AbstractPolicy<Vector, int>
{
protected:
	int n_actions;
	RandomNumberGenerator* rng;
public:
	RandomPolicy(int n_actions_, RandomNumberGenerator* rng_)
		: n_actions(n_actions_),
		  rng(rng_)
	{
	}
	virtual ~RandomPolicy()
	{
	}
	virtual int SelectAction()
	{
		return rng->discrete_uniform(n_actions);
	}
	virtual void Observe (const Vector& previous_state, const int& action, real r, const Vector& next_state) 
	{
	}
    virtual void Observe (real r, const Vector& next_state) 
	{
	}
	virtual void Reset() 
	{
	}
};
#endif
