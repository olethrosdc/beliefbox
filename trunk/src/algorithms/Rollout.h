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

#include "Environment.h"

/** A rollout.
 */
template<typename S, typename A, typename P>
class Rollout
{
public:
	S start_state; ///< initial state
	A start_action; ///< initial action
	S end_state; ///< current final state
	P* policy; ///< policy to be used for continuation
	Environment<S, A>* environment; ///< environment model
	real gamma; ///< discount factor
	int T; ///< total steps used
	real total_reward; ///< total reward received
	real discounted_reward; ///< discounted reward received
	bool running; ///< whether the rollout needs to be sampled again
	Rollout(const S& start_state_,
			const A& start_action_, 
			P* policy_, 
			Environment<S, A>* environment_,
			real gamma_) :
		start_state(start_state_),
		start_action(start_action_),
		policy(policy_),
		environment(environment_),
		gamma(gamma_)
	{
		end_state = start_state;
		total_reward = 0;
		discounted_reward = 0;
		running = false;
	}
	void Act(const A& a)
	{
		running = environment->Act(a);
		real reward = environment->getReward();
		//		printf("A: %d, r: %f ", a, reward);
		end_state = environment->getState();
		T++;
		total_reward += reward;
		discounted_reward += reward * gamma;
	}
	void Sample(const int period)
    {
        environment->Reset();
		environment->setState(start_state);
		running = true;
		for (int t=0; t<period; ++t) {
			if (!running) {
				break;
			}
			if (t==0) {
				Act(start_action);
			} else {
				Act(policy->SelectAction());
			}
		}
		printf("Length: %d, Total Reward: %f, State: ", T, total_reward);
		end_state.print(stdout);
	}
};


#endif
