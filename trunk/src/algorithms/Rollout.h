/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ROLLOUT_H
#define ROLLOUT_H

template<typename S, typename A, typename P>
class Rollout
{
public:
	S start_state;
	A action;
	S end_state;
	real gamma;
	int T;
	real total_reward;
	real discounted_reward;
	bool running;
	Rollout(S start_state_, A start_action_, P policy_, E environment_) :
		start_state(start_state_),
		start_action(start_action),
		policy(policy_),
		environment(environment_),
		gamma(gamma_)
	{
		end_state = start_state;
		total_reward = 0;
		discounted_reward = 0;
		running = false;
	}
	void Act(A a)
	{
		running = environment.Act(a);
		real reward = environment.getReward();
		end_state = environment.getState();
		T++;
		total_reward += reward;
		discounted_reward += reward * gamma;
	}
	void Sample(int period)
	{
		environment.setState(start_state);
		running = true;
		for (int t=0; t<period; ++t) {
			if (!running) {
				return;
			}
			Act(policy.SelectAction());
		}
	}
};


#endif
