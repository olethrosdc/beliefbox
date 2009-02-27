// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_MDP_H
#define RANDOM_MDP_H

#include "DiscreteMDP.h"
#include "Environment.h"
#include "RandomNumberGenerator.h"
#include <string>
#include <vector>



class RandomMDP : public DiscreteEnvironment
{
protected:
    uint n_states;
    uint n_actions;
    RandomNumberGenerator* rng;
public:
    RandomMDP(uint n_actions,
              uint n_states,
              real randomness,
              real step_value,
              real pit_value,
              real goal_value,
              RandomNumberGenerator* rng,
              bool termination=true);

    virtual ~RandomMDP();

    virtual DiscreteMDP* getMDP() const
    {
        return mdp;
    }

    /// put the environment in its natural state
    virtual void Reset();

    /// returns true if the action succeeds
    virtual bool Act(int action);

    /// Remove periods from MDP
    void AperiodicityTransform(real tau) 
    {
        mdp->AperiodicityTransform(tau);
    }
protected:
    DiscreteMDP* mdp;
    real** transitions;
    real* P_data;
    std::vector<Distribution*> rewards;
};

#endif
