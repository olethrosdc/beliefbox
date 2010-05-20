// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
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

POMDPGridworld::POMDPGridworld(RandomNumberGenerator* rng_,
                               const char* fname,
                               uint n_obs_,
                               real random_,
                               real pit_,
                               real goal_,
                               real step_)
    :  rng(rng_),
       observation(0),
       total_time(0),
       n_obs(n_obs_),
       random(random_), pit_value(pit_), goal_value(goal_), step_value(step_)
{
    n_actions = 4;


    // the number of obsrvations
    if (n_obs != 16 && n_obs != 2) {
        Serror ("Either 2 or 16 observations must be used\n");
    }

    CalculateDimensions(fname);

    // one state for the terminal state
    n_states = width * height + 1;
    terminal_state = n_states - 1;

    printf("# Making POMDPGridworld %d x %d from %s with %d states, %d obs, %f randomness, %f pit, %f goal, %f step rewards\n",
           height, width, fname, n_states, n_obs,
           random,
           pit_value,
           goal_value,
           step_value);
           
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
            case '.':
                grid[x][y] = GRID; 
                //n_states++;
                break;
            case '#':
                grid[x][y] = WALL; 
                break;
            case 'X':
                grid[x][y] = GOAL;
                //n_states++;
                break;
            case 'O':
                grid[x][y] = PIT;
                //n_states++;
                break;
            default: std::cerr << "Unknown maze element: '"
                               << line[x] << "'" << std::endl;
                exit(-1);
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
    
    
    pomdp = new DiscretePOMDP(n_states, n_obs, n_actions); // n_states?
    MakePOMDP();
    Reset();
}



POMDPGridworld::~POMDPGridworld() {
    for (uint i=0; i<rewards.size(); ++i) {
        delete rewards[i];
    }
    delete pomdp;
    //delete mdp;
}

void POMDPGridworld::CalculateDimensions(const char* fname)
{
    std::ifstream ifs(fname, std::ifstream::in);
    if (!ifs.is_open()) {
        Serror ("Could not open file %s", fname);
        exit(-1);
    }
    std::string line;
    height = 0;
    width = 0;
    while (getline(ifs, line)) {// && y<height) {
        if (!width) {
            width = line.length();
            if (!width) {
                Serror ("Empty first line\n");
                exit(-1);
            }
        } else if (line.length() != width) {
            Serror ("Line length (%ld) does not match width (%d)",
                    line.length(), width);
            exit(-1);
        }
        height++;
    }
    ifs.close();
}

void POMDPGridworld::MakePOMDP()
{
    //for (uint i=0; i
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
    } while(whatIs(x, y) != GRID || whatIs(x, y) == INVALID || whatIs(x, y) == WALL);
    //printf ("# Resetting to %d %d (%d)\n", x, y,  whatIs(x, y));
    ox = x;
    oy = y;
    switch(n_obs) {
    case 2:
        observation = 0;
        break;
    case 16:
        observation = CalculateObservation16obs();
        break;
    default:
        Serror("Invalid number of observations %d\n", n_obs);
    }
}

bool POMDPGridworld::Act(int action)
{
    observation = 0;

    if (state == terminal_state) {
        reward = 0;
        Reset();
        return false;
    }
    int x = ox;
    int y = oy;
    // Agent gets stuck with probability 'random'
    if (rng->uniform() >= random) {
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
    }
    reward = step_value;
    
    //if (rng->uniform() < random) {
    //reward = 1.0;
    //}
    //printf ("%d %d (%d) -> ", ox, oy, whatIs(ox, oy));
    //printf ("%d %d (%d) %d\n", x, y, whatIs(x, y), action);
    
    if (whatIs(x,y) != INVALID && whatIs(x,y) != WALL) {
        ox = x;
        oy = y;
        state = getStateFromCoords(x, y);
        if (n_obs == 2) {
            observation = 1;
        }
    } 

    if (whatIs(x,y) == GOAL || whatIs(x,y)== PIT) {
        state = terminal_state;
    }

    if (whatIs(x,y) == GOAL) {
        reward = goal_value;
    }

    if (whatIs(x,y) == PIT) {
        reward = pit_value;
    }
    //std::cout << "coords: " << x << ", " << y << std::endl;
    //if (state==terminal_state) {
        //std::cout << "t: " << total_time << " TERMINATE "
        //<< reward << std::endl;
    //return false;
        //}
    if (n_obs == 16) {
        observation = CalculateObservation16obs();
    }

    // Do not detect hit in some cases
    if (rng->uniform() < random) {
        observation = rng->discrete_uniform(n_obs);
    }
    return true;
}

void POMDPGridworld::Show()
{
    for (int y=0; y<(int) height; ++y) {
        for (int x=0; x<(int) width; ++x) {
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

int POMDPGridworld::CalculateObservation16obs()
{
    int obs = 0;

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
    assert(obs>=0 && obs<n_obs);
    return obs;
}
