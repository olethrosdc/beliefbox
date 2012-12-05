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
    uint ox, oy;
    int total_time;
    enum MapDirection {
        NORTH=0, SOUTH, EAST, WEST
    };
    enum MapElement {
        INVALID=-1, GRID, WALL, GOAL, PIT
    };
    DiscreteMDP* my_mdp;
    uint terminal_state;
    Gridworld(const char* fname,
              real random_ = 0.0,
              real pit_ = -1.0,
              real goal_ = 1.0,
              real step_ = -0.1);
    virtual ~Gridworld();

    static void GetMazeDimensions(const char* fname);
    
    virtual DiscreteMDP* getMDP() const;

    MapElement whatIs(int x, int y) const
    {
        if (x>=0 && y >=0 && x< (int) width && y < (int) height) {
            return grid[x][y];
        } else {
            return INVALID;
        }
    }
    virtual void Reset();
    virtual bool Act(const int& action);
    void Show();
    int getState(int x, int y) const
    {
        if (x>=0 && y >=0 && x< (int) width && y < (int) height) {
            return x + y*width;
        }
        return -1;
    };

    uint getWidth() const
    {
        return width;
    }

    uint getHeight() const
    {
        return height;
    }
    virtual const char* Name()
    {
        return "Gridworld";
    }
	virtual real getTransitionProbability(const int& state, const int& action, const int& next_state) const
	{
		return my_mdp->getTransitionProbability(state, action, next_state);
	}
    virtual real getExpectedReward(const int& state, const int& action) const
	{
		return my_mdp->getExpectedReward(state, action);
	}

protected:
    void CalculateDimensions(const char* fname);
    uint height;
    uint width;
    //uint n_aactions;
    real random;
    real pit_value;
    real goal_value;
    real step_value;
    std::vector< std::vector<MapElement> > grid;
    real** transitions;
    real* P_data;
    std::vector<Distribution*> rewards;
};

#endif
