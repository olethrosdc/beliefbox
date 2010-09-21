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

#ifndef RSAPI_H
#define RSAPI_H

#include "Rollout.h"
#include <vector>

template<typename R>
class StateSample
{
	std::vector<R> 
};

template<typename S, typename A>
class RSAPI
{
public:
	StateSet<S>
	A action;
	S end_state;
	real total_reward;
	real discounted_reward;
};


#endif
