/* -*- Mode: C++; -*- */
/* VER: $Id: Sampling.h,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MDP_DISTRIBUTION_H
#define MDP_DISTRIBUTION_H

#include "Object.h"
#include "MDP.h"
#include "Vector.h"
#include <vector>

/** A distribution over MDPs - these are observable MDPs, so the
	templatisation is over the observed variables: state and action.
	Rewards are always real numbers.
 */
template <typename StateType, typename ActionType>
class MDPDistribution : public Object
{
public:
	/// Destroy the MDP Distribution
	virtual ~MDPDistribution() {}

	/// Generate a sample from the distribution.
	virtual MDP<StateType, ActionType> generate() = 0;

	/// Register an observation, thus changing the distribution.
	virtual void observe (StateType s, ActionType a, real r, StateType s2) = 0;
};

/// Specialisation for distributions over discrete MDPs
typedef MDPDistribution<int, int> DiscreteMDPDistribution;

/// Specialisation for distributions over continuous state MDPs
typedef MDPDistribution<Vector, int> ContinuousStateMDPDistribution;

/// Specialisation for distributions over continuous action MDPs
typedef MDPDistribution<int, Vector> ContinuousActionMDPDistribution;

/// Specialisation for distributions over continuous MDPs
typedef MDPDistribution<Vector, Vector> ContinuousMDPDistribution;


#endif
