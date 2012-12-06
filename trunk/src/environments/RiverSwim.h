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

#ifndef RIVER_SWIM_H
#define RIVER_SWIM_H

#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "DiscreteMDP.h"

/// A simple environment to test exploration.
///
/// Adapted from Strehl and Littman 2006, "An analysis of model-based
/// interval estimation for Markov decision processes"
class RiverSwim : public DiscreteEnvironment
{
protected:
	real r_start; ///< reward at start state
	real r_end; ///< reward at end state
	real p_right; ///< probability of moving right, when trying to go right
	real p_left; ///< probability of moving left, when trying to go right
	real p_stuck; ///< probability of not moving, when trying to go right
public:
    DiscreteMDP* model; ///< internal model
    RiverSwim(int n = 6, real r_start_ = 0.0005, real r_end_ = 1.0);
    
    virtual ~RiverSwim();
    
    virtual void Reset();
    virtual bool Act(const int& action);

    virtual const char* Name()
    {
        return "River Swim";
    }

    virtual DiscreteMDP* getMDP() const;
};

#endif
