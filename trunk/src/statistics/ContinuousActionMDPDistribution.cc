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

#include "ContinuousActionMDPDistribution.h"

#if 0
ContinuousBanditDistribution::~ContinuousBanditDistribution()
{
}

ContinuousActionMDP ContinuousBanditDistribution::generate()
{
	return ContinuousActionMDP(transition_prior.generate(), reward_prior.generate());
}
void ContinuousBanditDistribution::observe (int s, Vector a, real r, int s2)
{
};

#endif
