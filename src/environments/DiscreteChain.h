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

/// A simple environment to test exploration.
///
/// This environment has been used in two different settings.
/// The One-Dimensional-Maze setting, and the Chain setting.
///
/// One-Dimensional-Maze
/// --------------------------------
/// Should be used with a discount of 0.75.
///
/// Returns a reward of 0 at every step,
/// and a reward of 1 at the end.
///
/// This is mainly useful when dealing with POMDPs.
///
/// Chain (default)
/// --------------
///
/// In this case, the environment is used to test exploration.
/// There is a slip probability of 0.2
/// The first action gives a reward of 0 apart from at the final state.
/// The second action always gives reward 0.2.
class DiscreteChain : public DiscreteEnvironment
{
protected:
	real slip, start, end;
public:
    DiscreteMDP* mdp;
    DiscreteChain(int n, real slip_ = 0.2, real start_ = 0.2, real end_ = 1.0);
    
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
