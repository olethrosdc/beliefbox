// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GRIDWORLD_H
#define GRIDWORLD_H

#include "DiscreteMDP.h"
#include <string>
#include <vector>

enum MapElement {
    GRID, WALL, GOAL, PIT
};

class Gridworld {
protected:
    uint height;
    uint width;
    uint n_actions;
    real random;
    real pit;
    std::vector< std::vector<MapElement> > grid;
    real** transitions;
    real* P_data;
    Distribution** rewards;
public:
    DiscreteMDP* mdp;
    Gridworld(char* fname,
	      uint height_,
	      uint width_,
	      uint n_actions_=4,
	      real random_=0.0,
	      real pit_=-100.0);
    
};

#endif
