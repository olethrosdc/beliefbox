// -*- Mode: c++ -*-
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDPCounts.h"

int main(void)
{
	int n_states = 5;
	int n_actions = 2;
	real dirichlet_mass = 0.5;
	enum DiscreteMDPCounts::RewardFamily reward_prior = DiscreteMDPCounts::BETA;
    DiscreteMDPCounts belief(n_states, n_actions, dirichlet_mass, reward_prior);

	DiscreteMDPCounts belief_copy(belief);

	belief.AddTransition(0, 1, 1, 1);
	belief_copy.AddTransition(0, 0, 0, 0);
	belief.ShowModelStatistics();
	belief_copy.ShowModelStatistics();

	DiscreteMDPCounts* clone = belief.Clone();

	delete clone;
	
	return 0;
}

