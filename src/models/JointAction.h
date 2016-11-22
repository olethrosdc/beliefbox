// -*- Mode: c++ -*-
// copyright (c) 2007-2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef JOINT_ACTION_H
#define JOINT_ACTION_H

#include "real.h"
#include "DiscreteStateSet.h"
#include "StateAction.h"
#include "HashCombine.h"
#include "debug.h"
#include <cstdio>
#include <map>
#include <unordered_map>


class DiscreteJointAction
{
public:
	int n_players; ///< number of players
	std::vector<int> n_actions; ///< number of actions per player
	std::vector<int> action; ///< joint action
	DiscreteJointAction(int n_players_,
						int n_actions_)
		: n_players(n_players_),
		  n_actions(n_players_),
		  action(n_players_)
	{
		for (int i=0; i<n_players; ++i) {
			action[i] = 0;
			n_actions[i] = n_actions_;
		}
	}
	
	/// Convert the action vector to an integer
	int to_int() const
	{
		int F = 1;
		int a = 0;
		for (int k=0; k<n_players; ++k) {
			a += action[k];
			F *= n_actions[k];
		}
		return a;
	}

	void set(const int& player, const int& action_)
	{
		assert(player >= 0 && player < n_players);
		assert(action >= 0 && action < n_actions);
		action[player] = action_;
	}
};


#endif
