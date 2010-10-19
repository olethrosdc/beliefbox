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

#ifndef CLASSIFIER_POLICY_H
#define CLASSIFIER_POLICY_H

#include "AbstractPolicy.h"

class ClassifierPolicy  : public AbstractPolicy<Vector, int>
{
protected:
	Classifier<Vector,int>* classifier;
public:
	ClassifierPolicy(Classifier<Vector,int>* classifier_);
	virtual ~ClassifierPolicy();
	virtual int SelectAction()
	{
		return classifier->Classify(state);
	}
	virtual void Observe (Vector& previous_state, ActionType& action, real r, StateType next_state) = 0;
    virtual void Observe (real r, StateType& next_state) = 0;
	virtual void Reset() = 0;
	virtual void SetState(StateType& state)
	{ 
		this->state = state;
	}
};
#endif
