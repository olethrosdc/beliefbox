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



class ContextBandit : public DiscreteEnvironment
{
protected:
    uint n_states;
    uint n_actions;
    uint terminal_state;
    RandomNumberGenerator* rng;
public:
    ContextBandit(uint n_actions,
                  uint n_features,
                  uint values_per_feature,
                  RandomNumberGenerator* rng);
    virtual ~ContextBandit()
    {
        for (uint i=0; i<rewards.size(); ++i) {
            delete rewards[i];
        }
        delete mdp;
    }
    /// put the environment in its natural state
    virtual void Reset();
    
    /// returns true if the action succeeds
    virtual bool Act(int action);

    virtual DiscreteMDP* getMDP() const
    {
        return mdp;
    }
    real getMean(int action);
    virtual const char* Name()
    {
        return "Context Bandit";
    }
protected:
    NormalDistribution normal;
    DiscreteMDP* mdp;
    std::vector<Distribution*> rewards;
};

#endif
