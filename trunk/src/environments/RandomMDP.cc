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

RandomMDP::RandomMDP(uint n_actions_,
                     uint n_states_,
                     real randomness,
                     real step_value,
                     real pit_value,
                     real goal_value,
                     bool termination) : n_states(n_states_),
                                         n_actions(n_actions_)
{
    
    // set up the mdp
    //std::cout << "Making the MPD\n";
    mdp = new DiscreteMDP (n_states, n_actions, NULL, NULL);
    
    //std::cout << "Setting up rewards\n";
    // set up rewards	
    
    SingularDistribution* step_reward = new SingularDistribution(step_value);
    SingularDistribution* pit_reward = new SingularDistribution(pit_value);
    SingularDistribution* zero_reward = new SingularDistribution(0.0);
    SingularDistribution* goal_reward = new SingularDistribution(goal_value);
    
    rewards.push_back(step_reward);
    rewards.push_back(pit_reward);
    rewards.push_back(zero_reward);
    rewards.push_back(goal_reward);
    
    //std::cout << "Assigning rewards\n";
    // first the terminal state rewards
    if (termination) {
        for (uint a=0; a<n_actions; ++a) {
            mdp->setRewardDistribution(n_states - 1, a, zero_reward);
        }
    }

    int terminal_state = n_states;
    if (termination) {
        terminal_state = n_states - 1;
    }
    for (uint s=0; s<terminal_state; ++s) { 
        for (uint a=0; a<n_actions; ++a) {
            mdp->setRewardDistribution(s, a, step_reward);
        }
    }

    int pit_state = (int) floor(((real) n_states - 1) * true_random(false));
    int goal_state = pit_state;
    while (goal_state == pit_state) {
        goal_state = (int) floor(((real) n_states - 1) * true_random(false));
    }



    for (uint a=0; a<n_actions; ++a) {
        mdp->setRewardDistribution(goal_state, a, goal_reward);
        mdp->setRewardDistribution(pit_state, a, pit_reward);
        if (termination) {
            mdp->setTransitionProbability (goal_state, a, terminal_state, 1.0);
            mdp->setTransitionProbability (pit_state, a, terminal_state, 1.0);
            mdp->setTransitionProbability (terminal_state, a, terminal_state, 1.0);
        }
    }


    //std::cout << "Assigning transition probabilities\n";
    // Step 1: set prior
    for (uint s=0; s<terminal_state; ++s) {   
        for (uint a=0; a<n_actions; ++a) {
            for (uint s2=0; s2<n_states - 1; s2++) {
                mdp->setTransitionProbability (s, a, s2, randomness);
            }
        }
    }
    
    // Step 2: add target states
    for (uint s=0; s<terminal_state; s++) {   
        for (uint a=0; a<n_actions; ++a) {
            int s2 = (int) floor(((real) n_states) * true_random(false));
            mdp->setTransitionProbability (s, a, s2, 1.0);
        }
    }

    // Step 3: normalize
    for (uint s=0; s<n_states; s++) {   
        for (uint a=0; a<n_actions; ++a) {
            real sum = 0.0;
            for (uint s2=0; s2<n_states; ++s2) {
                sum += mdp->getTransitionProbability (s, a, s2);
            }
            for (uint s2=0; s2<n_states; ++s2) {
                real p = mdp->getTransitionProbability (s, a, s2);
                mdp->setTransitionProbability (s, a, s2, p / sum);
            }
        }
    }
    fflush(stdout);
    mdp->Check();
}
	
RandomMDP::~RandomMDP() {
    for (uint i=0; i<rewards.size(); ++i) {
        delete rewards[i];
    }
    delete mdp;
}

/// put the environment in its natural state
void RandomMDP::Reset()
{
    state = (int) floor(urandom(0, n_states - 1));
    reward = 0;
}

/// returns true if the action succeeds, false if we are in a terminal state
bool RandomMDP::Act(int action)
{
    state = mdp->generateState(state, action);
    reward = mdp->generateReward(state, action);
    //reward = mdp->getExpectedReward(state, action);
    if (state != n_states - 1) {
        return true;
    }
    return false;  // environment has terminated
}
