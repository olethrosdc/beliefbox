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

#include "ContinuousActionMDP.h"

ContinuousBanditMDP::ContinuousBanditMDP(ContinuousActionTransitionDistribution& transition_distribution_, ContinuousActionRewardDistribution& reward_distribution_) : ContinuousActionMDP (transition_distribution_, reward_distribution_)
{

}

