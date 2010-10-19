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

#include "RSAPI.h"

RolloutState::RolloutState(Environment<Vector, int>& environment_,
						   Vector& start_state_)
	: environment(environment_),
	  start_state(start_state_)
{
	V_U = 0;
	V = 0;
}

RolloutState::~RolloutState()
{
	for (uint i=0; i<rollouts.size(); ++i) {
		delete rollouts[i];
	}
}

void RolloutState::NewRollout(AbstractPolicy<Vector, int> & policy, int action)
{
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout
		= new Rollout<Vector, int, AbstractPolicy<Vector, int> >(start_state, action, policy, environment, gamma);
	rollouts.push_back(rollout);
	
}	

void RolloutState::Sample(const int T)
{
	for (uint i=0; i<rollouts.size(); ++i) {
		rollouts[i]->Sample(T);
	}

}
