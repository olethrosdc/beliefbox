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

#ifndef DOUBLE_LOOP_H
#define DOUBLE_LOOP_H

#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "DiscreteMDP.h"

/// A simple environment to test exploration.
///
/// Adapted from Dearden, 1998, "Bayesian Q-Learning"
class DoubleLoop : public DiscreteEnvironment
{
protected:
	real r_left; ///< reward at start state
	real r_right; ///< reward at end state
public:
    DiscreteMDP* model; ///< internal model
    DoubleLoop(real r_left_ = 2.0, real r_right_ = 1.0);
    
    virtual ~DoubleLoop();
    
    virtual void Reset();
    virtual bool Act(int action);

    virtual const char* Name()
    {
        return "River Swim";
    }

    virtual DiscreteMDP* getMDP() const;
};

#endif
