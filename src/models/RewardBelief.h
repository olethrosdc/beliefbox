// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REWARD_BELIEF_H
#define REWARD_BELIEF_H

#include "RewardDistribution.h"
#include "Vector.h"
#include <vector>


/** A belief about rewards in discrete spaces.
    
 */
class DiscreteRewardBelief
{
protected:
    int n_states;
    int n_actions;
public:
    // only set up
    DiscreteRewardBelief(int n_states_, int n_actions_)
        : n_states(n_states_), n_actions(n_actions_)
    {}
        
    virtual ~DiscreteRewardBelief()
    {}

    virtual real Update(int state, int action, real reward) = 0;
};


/** A belief about a finite set of reward distributions in discrete
    spaces.
    
 */
class FiniteDiscreteRewardBelief
{
protected:
    std::vector<const DiscreteSpaceRewardDistribution*>& rewards; ///< reward distributions
    Vector probabilities; ///< vector of probabilites for each reward distribution
public:
    FiniteDiscreteRewardBelief(int n_states_, int n_actions_,
                               std::vector<const DiscreteSpaceRewardDistribution*>& rewards_,
                               Vector probabilities_)
        : DiscreteRewardBelief(n_states_, n_actions_),
          rewards(rewards_),
          probabilities(probabilites_)
    {
        assert(probabilites.Size() == reward
    }
        
    virtual ~FiniteDiscreteRewardBelief()
    {}
    virtual real Update(int state, int action, real reward);
};


