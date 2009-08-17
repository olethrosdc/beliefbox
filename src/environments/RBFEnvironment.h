// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RBF_ENVIRONMENT_H
#define ENVIRONMENT_H

#include "Environment.h"
#include "BasisSet.h"


class RBF_Environment : DiscreteEnvironment
{
protected:
    RBF_Environment<Vector, int>& environment;
public:
    RBF_Environment(Environment<Vector, int>& environment_) 
        : environment(environment_)
    {
    }
    virtual ~RBF_Environment() 
    {
    }

    /// put the environment in its natural state
    virtual void Reset() = 0;

    /// returns true if the action succeeds
    virtual bool Act(A action) = 0;

    virtual MDP<S, A>* getMDP() const
    {
        return NULL;
    }
    /// returns the current state
    S getState()
    {
        return state;
    }

    /// returns the current reward
    real getReward()
    {
        return reward;
    }

};

typedef Environment<int, int> DiscreteEnvironment;


#endif
