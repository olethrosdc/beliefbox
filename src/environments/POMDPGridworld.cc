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

#include "POMDPGridworld.h"
#include "Distribution.h"
#include <string>
#include <iostream>
#include <fstream>

POMDPGridworld::POMDPGridworld(const char* fname,
                     uint height_,
                     uint width_,
                     uint n_actions_,
                     real random_,
                     real pit_,
                     real goal_,
                     real step_)
    :  total_time(0),
       height(height_), width(width_), 
       random(random_), pit_value(pit_), goal_value(goal_), step_value(step_)
{
    n_states = width * height + 1; // plus a terminal state
    n_actions = n_actions_;
    terminal_state = n_states - 1;

    n_obs = 16;

    std::ifstream ifs(fname, std::ifstream::in);
    if (!ifs.is_open()) {
        Serror ("Could not open file %s", fname);
        exit(-1);
    }
    std::string line;
    grid.resize(width);
    for (uint i=0; i<width; i++) {
        grid[i].resize(height);
    }
    uint y = 0;
    while (getline(ifs, line)) {// && y<height) {
        if (line.length() != width) {
            Serror ("Line length (%ld) does not match width (%d)",
                    line.length(), width);
            exit(-1);
        }
        
        for (uint x=0; x<width; ++x) {
            switch (line[x]) {
            case '.': grid[x][y] = GRID; break;
            case '#': grid[x][y] = WALL; break;
            case 'X': grid[x][y] = GOAL; break;
            case 'O': grid[x][y] = PIT; break;
            default: std::cerr << "Unknown maze element\n"; exit(-1);
            }
        }
        y++;
    }				

    if (y < height) {
        std::cerr << "Only " << y << " lines read while accessing file "
                  << fname << std::endl;
        exit(-1);
    } else  if (y > height) {
        std::cerr << "Too many (" << y << ") lines read while accessing file "
                  << fname << std::endl;
        exit(-1);
    }
    
    Reset();
}

POMDPGridworld::~POMDPGridworld() {
    for (uint i=0; i<rewards.size(); ++i) {
        delete rewards[i];
    }
    //delete mdp;
}

void POMDPGridworld::Reset()
{
    total_time = 0;
    int x, y;
    int n_gridpoints = height*width;
    do {
        state = rand()%(n_gridpoints);
        x = state % height;
        y = (state - x) / width;
    } while(whatIs(x, y) != GRID);
    ox = x;
    oy = y;
}

bool POMDPGridworld::Act(int action)
{
    int x = ox;
    int y = oy;
    switch(action) {
    case NORTH:
        y--;
        break;
    case SOUTH:
        y++;
        break;
    case EAST:
        x++;
        break;
    case WEST:
        x--;
        break;
    }
    if (whatIs(x,y) != INVALID) {
        state = x;
    }
    if (state==terminal_state) {
        std::cout << "t: " << total_time << " TERMINATE "
                  << reward << std::endl;
        return false;
    }
    return true;
}

void POMDPGridworld::Show()
{
    for (int x=0; x<(int) width; ++x) {
        for (int y=0; y<(int) height; ++y) {
            if (x == ox && y == oy) {
                std::cout << "*";
                continue;
            }
            MapElement e = whatIs(x, y);
            switch (e) {
            case INVALID: std::cout << "!"; break;
            case GRID: std::cout << "."; break;
            case WALL: std::cout << "#"; break;
            case GOAL: std::cout << "X"; break;
            case PIT: std::cout << "O"; break;
            default: std::cout << "?"; break;
            }
        }
        std::cout << std::endl;
    }
}

int POMDPGridworld::getObservation()
{
    int obs = 0;
#if 0
    for (int y=oy-1; oy<=oy+1; y++) {
        for (int x = ox-1; x<=ox+1; x++) {
            if (y == oy && x == ox) 
                continue;
            MapElement e = whatIs(x, y);
            switch(e) {
            case GRID:
            case GOAL:
            case PIT:
                break;
            default:
                obs |= 1;
            }
            obs <<= 1;
        }
    }
#else
    if (whatIs(ox-1, oy)==WALL || whatIs(ox-1, oy)==INVALID) {
        obs |= 1;
    }

    if (whatIs(ox+1, oy)==WALL || whatIs(ox+1, oy)==INVALID) {
        obs |= 2;
    }

    if (whatIs(ox, oy-1)==WALL || whatIs(ox, oy-1)==INVALID) {
        obs |= 4;
    }

    if (whatIs(ox, oy+1)==WALL || whatIs(ox, oy+1)==INVALID) {
        obs |= 8;
    }
    
#endif
    assert(obs>=0 && obs<n_obs);
    return obs;
}
