// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_CHAIN_H
#define DISCRETE_CHAIN_H

#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "DiscreteMDP.h"

/// OneD Maze.
///
/// Should be used with a discount of 0.75.
///
/// Returns a reward of 0 at every step,
/// and a reward of 1 at the end.
class DiscreteChain : public DiscreteEnvironment
{
public:
    DiscreteMDP* mdp;
    DiscreteChain(int n);
    
    virtual ~DiscreteChain();
    
    virtual void Reset();
    virtual bool Act(int action);

    virtual const char* Name()
    {
        return "Discrete Chain";
    }

    virtual DiscreteMDP* getMDP() const;
};

#endif
