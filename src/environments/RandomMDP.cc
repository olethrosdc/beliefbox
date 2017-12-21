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

RandomMDP::RandomMDP(uint n_states_,
                     uint n_actions_,
                     real randomness_,
                     real step_value_,
                     real pit_value_,
                     real goal_value_,
                     RandomNumberGenerator* rng,
                     bool termination_)
    : DiscreteEnvironment(n_states_, n_actions_),
      randomness(randomness_),
      step_value(step_value_),
      pit_value(pit_value_),
      goal_value(goal_value_),
    termination(termination_)
{
    
    this->rng = rng;
    terminal_state = n_states;
    if (termination) {
        terminal_state = n_states - 1;
    }

    // set up the mdp
    std::cout << "Making the MPD\n";
    mdp = generateMDP();
    mdp->Check();
    Reset();
    printf ("state: %d/%d\n", state, n_states);
    printf("States: [%d, %d)\n", state_lower_bound, state_upper_bound);

}
	
DiscreteMDP* RandomMDP::generateMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP (n_states, n_actions, NULL);
    
    // first the terminal state rewards
    if (termination) {
        for (uint a=0; a<n_actions; ++a) {
            mdp-> setFixedReward(n_states - 1, a, 0.0);
        }
    }

    for (uint s=0; s<terminal_state; ++s) { 
        for (uint a=0; a<n_actions; ++a) {
            mdp->setFixedReward(s, a, step_value);
        } 
    }


    std::cout << "Setting pit/goal states\n";
    int pit_state = (int) floor(((real) n_states - 1) * rng->uniform());
    int goal_state = pit_state;
    while (goal_state == pit_state) {
        goal_state = (int) floor(((real) n_states - 1) * rng->uniform());
    }
    for (uint a=0; a<n_actions; ++a) {
        mdp->setFixedReward(goal_state, a, goal_value);
        mdp->setFixedReward(pit_state, a, pit_value);
        if (termination) {
            mdp->setTransitionProbability (goal_state, a, terminal_state, 1.0);
            mdp->setTransitionProbability (pit_state, a, terminal_state, 1.0);
            mdp->setTransitionProbability (terminal_state, a, terminal_state, 1.0);
        }
    }

    
    std::cout << "Assigning transition probabilities\n";
    
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
            int s2 = (int) floor(((real) n_states) * rng->uniform());
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

    return mdp;
}

RandomMDP::~RandomMDP() {
    delete mdp;
}

/// put the environment in its natural state
void RandomMDP::Reset()
{
    state = (int) rng->discrete_uniform(n_states);
    mdp->setState(state);
    reward = 0;
}

/// returns true if the action succeeds, false if we are in a terminal state
bool RandomMDP::Act(const int& action)
{
  reward = mdp->Act(action); //mdp->generateReward(state, action);
  state = mdp->getState();//mdp->generateState(state, action);
  if (state==terminal_state) {
      return false;
  } else {
    return true;  // we continue
  }
}
