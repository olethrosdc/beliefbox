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
	virtual void Observe (StateType previous_state, ActionType action, real r, StateType next_state) = 0;
    virtual void Observe (real r, StateType next_state) = 0;
	virtual void Reset() = 0;
	virtual void SetState(StateType state)
	{ 
		this->state = state;
	}
};


#endif
