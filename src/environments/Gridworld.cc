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

#include "Gridworld.h"
#include "Distribution.h"
#include <string>
#include <iostream>
#include <fstream>

Gridworld::Gridworld(const char* fname,
                     uint height_,
                     uint width_,
                     uint n_actions_,
                     real random_,
                     real pit_,
                     real goal_,
                     real step_)
    :  height(height_), width(width_), n_actions(n_actions_),
       random(random_), pit_value(pit_), goal_value(goal_), step_value(step_)
{
    uint n_states = width * height + 1; // plus a terminal state

    std::ifstream ifs(fname, std::ifstream::in);
    std::string line;
    grid.resize(width);
    for (uint i=0; i<width; i++) {
        grid[i].resize(height);
    }
    uint y = 0;
    while (getline(ifs, line) && y<height) {
        assert (line.length() == width);
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


    // set up the mdp

    mdp = new DiscreteMDP (n_states, n_actions, NULL, NULL);

    // set up rewards		
    SingularDistribution step_reward(step_value);
    SingularDistribution pit_reward(pit_value);
    SingularDistribution zero_reward(0.0);
    SingularDistribution goal_reward(goal_value);
    
    // first the terminal state rewards
    for (uint a=0; a<n_actions; ++a) {
        mdp->setRewardDistribution(n_states - 1, a, &zero_reward);
    }
    // then the others.
    for (uint x=0; x<width; ++x) {
        for (uint y=0; y<height; ++y) {
            for (uint a=0; a<n_actions; ++a) {
                int s = getState(x,y);
                switch(whatIs(x,y)) {
                case GRID:
                    mdp->setRewardDistribution(s, a, &step_reward);
                    break;
                case WALL:
                    mdp->setRewardDistribution(s, a, &zero_reward);
                    break;
                case GOAL:
                    mdp->setRewardDistribution(s, a, &goal_reward);
                    break;
                case PIT:
                    mdp->setRewardDistribution(s, a, &pit_reward);
                    break;
                default:
                    std::cerr << "Unknown grid point type\n";
                    exit(-1);
                }
            }
        }
    }
    
    // set up transitions
    // first the terminal state
    for (uint a=0; a<n_actions; ++a) {
        mdp->setTransitionProbability (n_states-1, a, n_states-1, 1.0);
        for (uint s=0; s<n_states - 1; s++) {
            mdp->setTransitionProbability (n_states-1, a, s, 0.0);
            mdp->setTransitionProbability (s, a, n_states-1, 0.0);
        }
    }

    // then all the other states
    // Step 1: clear
    for (uint s=0; s<n_states -1; s++) {   
        for (uint a=0; a<n_actions; ++a) {
            for (uint s2=0; s2<n_states - 1; s2++) {
                mdp->setTransitionProbability (s, a, s2, 0.0);
            }
        }
    }

    // Step 2: fill
    for (uint x=0; x<width; ++x) {
        for (uint y=0; y<height; ++y) {
            uint s = getState(x, y);
            MapElement element = whatIs(x, y);
            if (element == WALL) {
                for (uint a=0; a<n_actions; ++a) {
                    mdp->setTransitionProbability (s, a, s, 1.0);
                }
                continue;
            } else if (element == GOAL || element == PIT) {
                std::cout << "TERMINATE: " << s << std::endl;
                for (uint a=0; a<n_actions; ++a) {
                    mdp->setTransitionProbability (s, a, n_states-1, 1.0);
                }
                continue;
            } else if (element == INVALID) {
                std::cerr << "Invalid element\n";
                exit(-1);
            }

            int num = 4;
            // the hardest part is checking walls
            bool Nd = true;
            bool Sd = true;
            bool Wd = true;
            bool Ed = true;
            int Es = getState(x + 1, y);
            int Ws = getState(x - 1, y);
            int Ns = getState(x, y - 1);
            int Ss = getState(x, y + 1);

            if (x==0 || whatIs(x-1, y) == WALL)  {
                Ws = s;
                Wd = false;
                num--;
            } 
            if (x==width-1 || whatIs(x+1, y) == WALL)  {
                Es = s;
                Ed = false;
                num--;
            } 
            if (y==0 || whatIs(x, y-1) == WALL)  {
                Ns = s;
                Nd = false;
                num--;
            } 
            if (y==height-1 || whatIs(x, y+1) == WALL)  {
                Ss = s;
                Sd = false;
                num--;
            } 

            real theta = random / (real) num;
            for (uint a=0; a<n_actions; ++a) {
                if (Ed) {
                    mdp->setTransitionProbability (s, a, Es, theta);
                }
                if (Wd) {
                    mdp->setTransitionProbability (s, a, Ws, theta);
                }
                if (Nd) {
                    mdp->setTransitionProbability (s, a, Ns, theta);
                }
                if (Sd) {
                    mdp->setTransitionProbability (s, a, Ss, theta);
                }
                switch(a) {
                case NORTH:
                    mdp->setTransitionProbability (s, a, Ns, 1 - random + theta);
                    break;
                case SOUTH:
                    mdp->setTransitionProbability (s, a, Ss, 1 - random + theta);
                    break;
                case EAST:
                    mdp->setTransitionProbability (s, a, Es, 1 - random + theta);
                    break;
                case WEST:
                    mdp->setTransitionProbability (s, a, Ws, 1 - random + theta);
                    break;
                }
            }
        }
    }
    //mdp->Check();
}

void Gridworld::Reset()
{
    int x, y;
    int n_gridpoints = height*width;
    do {
        state = rand()%(n_gridpoints);
        x = state % height;
        y = state - x;
    } while(whatIs(x, y) != GRID);
}

bool Gridworld::Act(int a)
{
    state = mdp->generateState(state, a);
    return true;
}

void Gridworld::Show()
{
    for (uint x=0; x<width; ++x) {
        for (uint y=0; y<height; ++y) {
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
