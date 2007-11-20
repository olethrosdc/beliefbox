/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTINUOUS_BANDIT_POLICY_H
#define CONTINUOUS_BANDIT_POLICY_H

class AbstractContinuousBanditPolicy : public AbstractContinuousActionPolicy
{
	virtual ~AbstractContinuousBanditPolicy();
};

class IntervalSamplingContinuousBanditPolicy : public AbstractContinuousBanditPolicy
{
	virtual ~IntervalSamplingContinuousBanditPolicy();
	virtual Vector SelectAction();
	virtual void Observe (int previous_state,
						  Vector& action,
						  real r,
						  int next_state);
	virtual void Reset();
};

#endif
