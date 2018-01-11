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
                             RandomNumberGenerator* rng_,
                             bool normal)
    : DiscreteEnvironment(n_states_, n_actions_),
      rng(rng_)
{ 
    logmsg ("Making bandit with %d states, %d actions\n", n_states, n_actions);
    mdp = new DiscreteMDP (n_states, n_actions, NULL);
    assert(rng);
    state = 0;
    reward = 0;
    // setup rewards
    for (uint s=0; s<n_states; ++s) {
        for (uint a=0; a<n_actions; a++) {
            Distribution* reward_dist;
            if (normal) {
                reward_dist = new NormalDistribution(normal_prior.generate(), 1.0);
            } else {
                reward_dist = new BernoulliDistribution(urandom());
            }
            mdp->setRewardDistribution(s, a, reward_dist);
            rewards.push_back(reward_dist);
            printf ("%d %d %f # E[r|s,a]\n",
                    s, a,
                    mdp->getExpectedReward(s, a));
        }
    }

    // set up transitions uniformly
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
    Reset();
}


/// put the environment in its natural state
void ContextBandit::Reset()
{
    state = (int) rng->discrete_uniform(n_states);
    reward = 0;
    mdp->setState(state);
}


/// returns true if the action succeeds, false if we are in a terminal state
bool ContextBandit::Act(const int& action)
{
    bool action_ok = mdp->Act(action);
    state = mdp->getState();
    reward = mdp->getReward();
    return action_ok;  // we continue
}


ContextBandit::~ContextBandit()
{
    for (uint i=0; i<rewards.size(); ++i) {
        delete rewards[i];
    }
    delete mdp;
}

DiscreteMDP* ContextBandit::getMDP() const
{
    return new DiscreteMDP(*mdp);
}

const char* ContextBandit::Name()
{
    return "Context Bandit";
}


