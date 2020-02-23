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
#include "SingularDistribution.h"
#include <string>
#include <iostream>
#include <fstream>

Gridworld::Gridworld(const char* fname,
                     real random_,
                     real pit_,
                     real goal_,
                     real step_)
    : total_time(0), start_state(false),
    random(random_), pit_value(pit_), goal_value(goal_), step_value(step_)
{
    n_actions = 4;
    
    CalculateDimensions(fname);
    n_states = width * height + 1; // plus a terminal state
    state_upper_bound = n_states;
    terminal_state = n_states - 1;
    //reward = 0.0;
    logmsg("Gridworld, states: %d, random: %f, pit: %f, goal: %f, step: %f\n",
           n_states,
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
                    (long int) line.length(),  width);
            exit(-1);
        }
        
        for (uint x=0; x<width; ++x) {
            switch (line[x]) {
            case '.': grid[x][y] = GRID; break;
            case '#': grid[x][y] = WALL; break;
            case 'X': grid[x][y] = GOAL; break;
            case 'O': grid[x][y] = PIT; break;
			case 'S': grid[x][y] = START; start_state_x = x; start_state_y = y; start_state = true; break;
            default: std::cerr << "Unknown maze element (" << line[x] << ")\n"; exit(-1);
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
    
    // set up the mdp
    
    my_mdp = getMDP();
    Reset();
}

DiscreteMDP* Gridworld::getMDP()  const
{
    DiscreteMDP* mdp = new DiscreteMDP (n_states, n_actions, NULL);

    // set up rewards		
#if 0
    SingularDistribution* step_reward = new SingularDistribution(step_value);
    SingularDistribution* pit_reward = new SingularDistribution(pit_value);
    SingularDistribution* zero_reward = new SingularDistribution(0.0);
    SingularDistribution* goal_reward = new SingularDistribution(goal_value);
    
#if 0
    std::cout << "step:" << step_value
              << " pit:" << pit_value
              << " goal:" << goal_value
              << std::endl;
#endif
    rewards.push_back(step_reward);
    rewards.push_back(pit_reward);
    rewards.push_back(zero_reward);
    rewards.push_back(goal_reward);
#endif

    // first the terminal state rewards
    for (uint a=0; a<n_actions; ++a) {
        mdp->setFixedReward(terminal_state, a, 0.0);
    }
    // then the others.
    for (uint x=0; x<width; ++x) {
        for (uint y=0; y<height; ++y) {
            for (uint a=0; a<n_actions; ++a) {
                int s = getState(x,y);
                switch(whatIs(x,y)) {
                case GRID:
                    mdp->setFixedReward(s, a, step_value);
                    break;
                case WALL:
                    mdp->setFixedReward(s, a, 0.0);
                    break;
                case GOAL:
                    mdp->setFixedReward(s, a, goal_value);
                    break;
                case PIT:
                    mdp->setFixedReward(s, a, pit_value);
                    break;
				case START:
                    mdp->setFixedReward(s, a, 0);
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
        mdp->setTransitionProbability (terminal_state, a, terminal_state, 1.0);
        for (uint s=0; s<terminal_state; s++) {
            mdp->setTransitionProbability (terminal_state, a, s, 0.0);
            //mdp->setTransitionProbability (s, a, terminal_state, 0.0);
        }
    }


    // then all the other states
#if 0
    // Step 1: initialise them all to be self-looping
    for (uint s=0; s<n_states -1; s++) {   
        for (uint a=0; a<n_actions; ++a) {
            mdp->setTransitionProbability (s, a, s, 1.0);
        }
    }
#endif

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
                //std::cout << "TERMINATE: " << s << " " << element << std::endl;
                for (uint a=0; a<n_actions; ++a) {
                    mdp->setTransitionProbability (s, a, terminal_state, 1.0);
                }
                continue;
            } else if (element == INVALID) {
                std::cerr << "Invalid element\n";
                exit(-1);
            } 

            int num = 4;

            // the hardest part is checking walls
            // check whether movement in a particular direction is possible
            bool Nd = true; 
            bool Sd = true;
            bool Wd = true;
            bool Ed = true;

            // Collect the types of states of all different directions
            int Es = getState(x + 1, y);
            int Ws = getState(x - 1, y);
            int Ns = getState(x, y - 1);
            int Ss = getState(x, y + 1);

            // collect movements
            if (x<=0 || whatIs(x-1, y) == WALL || whatIs(x-1, y) == INVALID)  {
                Wd = false;
                num--;
            } 
            if (x>=width-1 || whatIs(x+1, y) == WALL || whatIs(x+1, y) == INVALID)  {
                Ed = false;
                num--;
            } 
            if (y<=0 || whatIs(x, y-1) == WALL || whatIs(x, y-1) == INVALID)  {
                Nd = false;
                num--;
            } 
            if (y>=height-1 || whatIs(x, y+1) == WALL || whatIs(x, y+1) == WALL)  {
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
                    if (Nd) {
                        mdp->setTransitionProbability (s, a, Ns, 1 - random + theta);
                        //printf("a: %d: %d -> %d\n", a, s, Ns);
                    } else {
                        mdp->setTransitionProbability (s, a, s, 1 - random + theta);
                    }
                    break;
                case SOUTH:
                    if (Sd) {
                        mdp->setTransitionProbability (s, a, Ss, 1 - random + theta);
                        //printf("a: %d: %d -> %d\n", a, s, Ss);
                    } else {
                        mdp->setTransitionProbability (s, a, s, 1 - random + theta);
                    }
                    break;
                case EAST:
                    if (Ed) {
                        mdp->setTransitionProbability (s, a, Es, 1 - random + theta);
                        //printf("a: %d: %d -> %d\n", a, s, Es);
                    } else {
                        mdp->setTransitionProbability (s, a, s, 1 - random + theta);
                    }
                    break;
                case WEST:
                    if (Wd) {
                        mdp->setTransitionProbability (s, a, Ws, 1 - random + theta);
                        //printf("a: %d: %d -> %d\n", a, s, Ws);
                    } else {
                        mdp->setTransitionProbability (s, a, s, 1 - random + theta);
                    }
                    break;
                }
            }
        }
    }

    for (uint s=0; s<n_states; ++s) {
        for (uint a=0; a<n_actions; ++a) {
            real sum = 0;
            for (uint s2=0; s2<n_states; ++s2) {
                sum += mdp->getTransitionProbability (s, a, s2);
            }
            //printf ("sum: %f -> ", sum);
            real isum = 1.0 / sum;
            sum = 0;
            for (uint s2=0; s2<n_states; ++s2) {
                real p = mdp->getTransitionProbability (s, a, s2);
                if (p > 0) {
                    mdp->setTransitionProbability (s, a, s2, p * isum);
                    sum += mdp->getTransitionProbability (s, a, s2);
                }
            }
            //printf (" %f\n", sum);
        }
    }

    //mdp->ShowModel();
    mdp->Check();
    return mdp;
}

Gridworld::~Gridworld() {
    for (uint i=0; i<rewards.size(); ++i) {
        delete rewards[i];
    }
    delete my_mdp;
}

void Gridworld::Reset()
{
    total_time = 0;
    int x, y;
    int n_gridpoints = height*width;
	if (start_state) {
		
		x = start_state_x;
		y = start_state_y;
		state = getState(x, y);
	} else {
		do {
        state = rand()%(n_gridpoints);
        x = state % height;
        y = (state - x) / width;
		} while(whatIs(x, y) != GRID);
	}
    ox = x;
    oy = y;
    reward = 0.0;
    my_mdp->Reset(state);
}

bool Gridworld::Act(const int& action)
{
    int x = state % width;
    int y = (state - x) / width;
    assert(state == getState(x, y));
    assert(whatIs(x,y) != INVALID && whatIs(x, y) != WALL);
    //std::cout << "(" << x << ", "<< y << ")" << " [" << whatIs(x,y) << "] a: " << action;
    total_time++;
    //real prev_reward = reward;
    reward = my_mdp->generateReward(state, action);
    state = my_mdp->generateState(state, action);
    x = state % width;
    y = (state - x) / width;

    assert((state == (int) terminal_state) || (whatIs(x,y) != INVALID && whatIs(x, y) != WALL));

    ox = x;
    oy = y;

	//	std::cout << " -> (" << x << ", "<< y << ")" << " [" << whatIs(x,y)
	//<< "] s: " << state << "r: " << reward << std::endl;
    //Show();

    if (state==(int) terminal_state) {
        return false;
    }
    return true;
}

void Gridworld::Show()
{


    for (uint y=0; y<height; ++y) {
        for (uint x=0; x<width; ++x) {
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
			case START: std::cout << "S"; break;
            default: std::cout << "?"; break;
            }
        }
        std::cout << std::endl;
    }
}

void Gridworld::CalculateDimensions(const char* fname)
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
