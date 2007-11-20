/* -*- Mode: C++; -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DBNMDP_H
#define DBNMDP_H


#include "MDPDistribution.h"

/** A distribution over MDPs - these are observable MDPs, so the
	templatisation is over the observed variables: state and action.
	Rewards are always real numbers.
 */
class DBNMDP :  public DiscreteMDPDistribution
{
public:
	/// Destroy the MDP Distribution
	virtual ~DBNMDP() {}

	/// Generate a sample from the distribution.
	virtual MDP<int, int> generate();

	/// Register an observation, thus changing the distribution.
	virtual void observe (int s, int a, real r, int s2) = 0;
};


#endif
