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
#include "Environment.h"
#include <string>
#include <vector>



class Gridworld : public Environment<int, int> {
public:
    enum MapDirection {
        NORTH=0, SOUTH, EAST, WEST
    };
    enum MapElement {
        INVALID=-1, GRID, WALL, GOAL, PIT
    };
    DiscreteMDP* mdp;

    Gridworld(const char* fname,
              uint height_,
              uint width_,
              uint n_actions_ = 4,
              real random_ = 0.0,
              real pit_ = -100.0,
              real goal_ = 0.0,
              real step_ = -1.0);
    virtual ~Gridworld()
    {
    }

    virtual DiscreteMDP* getMDP() const
    {
        return mdp;
    }

    MapElement whatIs(int x, int y)
    {
        if (x>=0 && y >=0 && x< (int) width && y < (int) height) {
            return grid[x][y];
        } else {
            return INVALID;
        }
    }
    virtual void Reset();
    virtual bool Act(int a);
    void Show();
    int getState(int x, int y)
    {
        if (x>=0 && y >=0 && x< (int) width && y < (int) height) {
            return x + y*width;
        }
        return -1;
    };
protected:
    uint height;
    uint width;
    uint n_actions;
    real random;
    real pit_value;
    real goal_value;
    real step_value;
    std::vector< std::vector<MapElement> > grid;
    real** transitions;
    real* P_data;
    Distribution** rewards;
};

#endif
