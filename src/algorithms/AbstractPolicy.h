/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ABSTRACT_POLICY_H
#define ABSTRACT_POLICY_H

#include <vector>
#include "real.h"
#include "MDPDistribution.h"

template <typename StateType, typename ActionType>
class AbstractPolicy
{
public:	
	StateType state;
	virtual ~AbstractPolicy(){};
	virtual ActionType SelectAction() = 0;
	virtual void Observe (const StateType& previous_state, const ActionType& action, real r, const StateType& next_state) = 0;
    virtual void Observe (real r, const StateType& next_state) = 0;
	virtual void Reset() = 0;
	virtual void setState(const StateType& state)
	{ 
		this->state = state;
	}
    virtual void setEpsilonGreedy(real epsilon_)
    {
        Serror("Not implented\n");
        exit(-1);
    }
	/// Update with gradient
	///
	/// TODO: Refactor PolicyGrdient.cc to use this function?
	virtual void GradientUpdate()
	{
	}
	
};


#endif
