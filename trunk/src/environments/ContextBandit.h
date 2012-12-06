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

#ifndef CONTEXT_BANDIT_H
#define CONTEXT_BANDIT_H

#include "DiscreteMDP.h"
#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "NormalDistribution.h"
#include <string>
#include <vector>


/** Context bandit
    
    In this MDP, the probability of the next state does not depend on
    the previous or the previous action. It is in fact drawn from a
    stochastically generated multinomial distribution.
    
    The rewards are randomly initialised, so that each state-action
    combination has a reward distribution of either Normal(-1,1), or
    Normal(1,1), but it is not known which.

 */
class ContextBandit : public DiscreteEnvironment
{
protected:
    RandomNumberGenerator* rng;
public:
    ContextBandit(uint n_states_,
                  uint n_actions_,
                  RandomNumberGenerator* rng_,
                  bool normal = true);
    virtual ~ContextBandit();

    /// put the environment in its natural state
    virtual void Reset();
    
    /// returns true if the action succeeds
    virtual bool Act(const int& action);

    /// Get the MDP
    virtual DiscreteMDP* getMDP() const;

    /// Get the name of the environment
    virtual const char* Name();

protected:
    NormalDistribution normal_prior;
    DiscreteMDP* mdp;
    std::vector<Distribution*> rewards;
};

#endif
