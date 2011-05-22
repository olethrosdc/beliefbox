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

#ifndef REWARD_POLICY_BELIEF_H
#define REWARD_POLICY_BELIEF_H

#include "PolicyBelief.h"


/** Prior on epsilon optimality of policies */
class RewardPolicyBelief
{
public:
    DirichletPolicyBelief policy_belief;
    const DiscreteSpaceRewardDistribution& rewards;
    RewardPolicyBelief(int n_states, int n_actions,
                       DiscreteSpaceRewardDistribution& rewards)
        : DirichletPolicyBelief(n_states, n_actions),
          
    {
        
    }
    
};

#endif
