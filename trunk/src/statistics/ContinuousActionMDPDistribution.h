/* -*- Mode: C++; -*- */
/* $Revision$ */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTINUOUS_ACTION_MDP_DISTRIBUTION_H
#define CONTINUOUS_ACTION_MDP_DISTRIBUTION_H

#include "ContinuousActionMDP.h"
#include "MDPDistribution.h"
#include "Vector.h"

/// Specialisation for distributions over continuous action MDPs
typedef MDPDistribution<int, Vector> ContinuousActionMDPDistribution;

/* Continuous bandit MDP distribution.
 *
 * This has only one state.
 */
class ContinuousBanditDistribution : public ContinuousActionMDPDistribution
{
protected:
	
public:
	virtual ~ContinuousBanditDistribution();
	virtual ContinuousActionMDP generate();
	virtual void observe (int s, Vector a, real r, int s2) = 0;
};

#endif
