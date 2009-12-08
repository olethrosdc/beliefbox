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

#include "DiscreteEnvironment.h"
#include "RandomNumberGenerator.h"

/// OneD Maze
///
/// Should be used with a discount of 0.75
class OneDMaze : public DiscreteEnvironment
{
protected:
    int hidden_state;
    int n_hidden_states;
    RandomNumberGenerator* rng;

    OneDMaze();

    virtual ~OneDMaze();
    
    virtual void Reset();
    virtual bool Act(A action);

    virtual const char* Name()
    {
        return "OneD Maze";
    }

};

#endif
