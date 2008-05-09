// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "RandomMDP.h"
#include "Distribution.h"
#include "Random.h"
#include <string>
#include <iostream>
#include <fstream>

RandomMDP::RandomMDP(uint n_actions,
                     uint n_states,
                     real randomness,
		     real step_value,
		     real pit_value,
		     real goal_value)
{
    
    // set up the mdp
    mdp = new DiscreteMDP (n_states, n_actions, NULL, NULL);

    // set up rewards		
    SingularDistribution step_reward(step_value);
    SingularDistribution pit_reward(pit_value);
    SingularDistribution zero_reward(0.0);
    SingularDistribution goal_reward(goal_value);
    
    // first the terminal state rewards
    for (uint a=0; a<n_actions; ++a) {
        mdp->setRewardDistribution(n_states - 1, a, &zero_reward);
    }

    for (uint s=0; s<n_states; s++) { 
	for (uint a=0; a<n_actions; ++a) {
	    mdp->setRewardDistribution(s, a, &step_reward);
	}
    }
    int pit_state = (int) floor(((real) n_states) * true_random());
    int goal_state = (int) floor(((real) n_states) * true_random());

    for (uint a=0; a<n_actions; ++a) {
	mdp->setRewardDistribution(pit_state, a, &pit_reward);
	mdp->setRewardDistribution(goal_state, a, &goal_reward);
    }

    // Step 1: clear
    for (uint s=0; s<n_states -1; s++) {   
        for (uint a=0; a<n_actions; ++a) {
            for (uint s2=0; s2<n_states - 1; s2++) {
                mdp->setTransitionProbability (s, a, s2, randomness);
            }
        }
    }
    
    // Step 2: add target states
    for (uint s=0; s<n_states -1; s++) {   
        for (uint a=0; a<n_actions; ++a) {
	    int s2 = (int) floor(((real) n_states) * true_random());
	    mdp->setTransitionProbability (s, a, s2, 1.0);
        }
    }

    // Step 3: add prior and normalize
    for (uint s=0; s<n_states -1; s++) {   
        for (uint a=0; a<n_actions; ++a) {
	    int s2 = (int) floor(((real) n_states) * true_random());
	    mdp->setTransitionProbability (s, a, s2, 1.0);
        }
    }


    mdp->Check();
}


