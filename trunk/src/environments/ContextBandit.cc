// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ContextBandit.h"


ContextBandit::ContextBandit(uint n_actions_,
                             uint n_features,
                             unit values_per_feature,
                             RandomNumberGenerator* rng): actions(n_actions_)
{ 
    this->rng = rng;
    // this bandit is binary
    n_states = (int) round(pow(values_per_feature, n_features));
}

/// put the environment in its natural state
void RandomMDP::Reset()
{
    state = (int) rng->discrete_uniform(n_states);
    reward = 0;
}

/// returns true if the action succeeds, false if we are in a terminal state
bool RandomMDP::Act(int action)
{
    int f0 = state & 7;
    int f1 = (state >> 3) & 7;
    int f2 = (state >> 6) & 7;
   	real mean;
	real sigma = 2;

	if ((f1 % 4) == 1)
	{
		if (f2 % 2)
		{
			mean = 1;
		}
		else
		{
			mean = -1;
		}
	}
	else
	{
		if (f0 % 2)
		{
			mean = 1;
		}
		else
		{
			mean = -1;
		}
	}
	if (action)   // a = 0 or a = 1
	{
		mean = -mean;
	}

    state = (int) rng->discrete_uniform(n_states);
    normal.setVariance(sigma*sigma);
    normal.setMean(mean);
    reward = normal.generate();
    //reward = mdp->getExpectedReward(state, action);
    if (state == terminal_state) {
        return false; // terminated
    }
    return true;  // we continue
}
