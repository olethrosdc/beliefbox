// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef OneD_MAZE_H
#define OneD_MAZE_H

#include "Environment.h"
#include "RandomNumberGenerator.h"

/// OneD Maze.
///
/// Should be used with a discount of 0.75.
///
/// Returns a reward of 0 at every step,
/// and a reward of 1 at the end.
class OneDMaze : public DiscreteEnvironment
{
protected:
    int hidden_state;
    int n_hidden_states;
    RandomNumberGenerator* rng;
public:
    OneDMaze(int n_hidden_states_, RandomNumberGenerator* rng_);

    virtual ~OneDMaze();
    
    virtual void Reset();
    virtual bool Act(const int& action);

    virtual const char* Name()
    {
        return "OneD Maze";
    }

};

#endif
