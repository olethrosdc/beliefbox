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

#ifndef RSAPI_H
#define RSAPI_H

#include "StateSet.h"
#include "Rollout.h"
#include "Vector.h"
#include <vector>

class RolloutState
{
public:
	Vector start_state;
	real gamma;
	std::vector<Rollout<Vector, int, AbstractPolicy<Vector, int> >* > rollouts;
};

template <typename S, typename A, typename P>
class RSAPI
{
public:
	std::vector<RolloutState*> states;
};

class void_RSAPI
{
public:
	std::vector<void*> states;
	
};

#endif
