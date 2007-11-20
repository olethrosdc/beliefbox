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

#include "MDP.h"
#include "Vector.h"

typedef MDP<int, Vector> ContinuousActionMDP;

typedef TransitionDistribution<int, Vector> ContinuousActionTransitionDistribution;
typedef RewardDistribution<int, Vector> ContinuousActionRewardDistribution;


/** Continuous bandit MDP.
 *
 * In this type of MDP there are bandits.
 * Hooray.
 * 
 */
class ContinuousBanditMDP : public ContinuousActionMDP
{

public:
	ContinuousBanditMDP(ContinuousActionTransitionDistribution& transition_distribution_, ContinuousActionRewardDistribution& reward_distribution_);
	virtual ~ContinuousBanditMDP();
};
