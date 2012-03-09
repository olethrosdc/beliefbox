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


ContextBandit::ContextBandit(uint n_states_,
                             uint n_actions_,
                             RandomNumberGenerator* rng_)
    : DiscreteEnvironment(n_states_, n_actions_),
      rng(rng_)
{ 
    printf ("Making bandit!\n");
    mdp = new DiscreteMDP (n_states, n_actions, NULL);
    assert(rng);
    state = 0;
    reward = 0;
    // setup rewards
    for (uint s=0; s<n_states; ++s) {
        for (uint a=0; a<n_actions; a++) {
            NormalDistribution* reward_dist
                = new NormalDistribution(urandom(), 1.0);
            mdp->setRewardDistribution(s, a, reward_dist);
            //rewards.push_back(reward_dist);
        }
    }

    // set up transitions
    real Pr = 1.0 / (real) n_states;
    for (uint s=0; s<n_states; s++) {
        for (uint a=0; a<n_actions; a++) {
            for (uint s2=0; s2<n_states; s2++) {
                //printf ("%d %d %d / %d", s, a, s2, n_states);
                mdp->setTransitionProbability(s, a, s2, Pr);
            }
        }
    }
    mdp->Check();
}


/// put the environment in its natural state
void ContextBandit::Reset()
{
    state = (int) rng->discrete_uniform(n_states);
    reward = 0;
}


/// returns true if the action succeeds, false if we are in a terminal state
bool ContextBandit::Act(int action)
{
    reward = mdp->Act(action);
    return true;  // we continue
}
