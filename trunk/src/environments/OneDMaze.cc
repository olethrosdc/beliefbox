/ -*- Mode: c++ -*-
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

#include "OneDMaze.h"
#include "debug.h"

OneDMaze::OneDMaze(int n_hidden_states_, RandomNumberGenerator* rng_) : Environment(2, 2), n_hidden_states(n_hidden_states_), rng(rng_)
{
    Reset();
}

OneDMaze::~OneDMaze()
{
}

void OneDMaze::Reset()
{
    hidden_state = 1 + rng->discrete_uniform(n_hidden_states - 1);
    reward = 0;
    state = 0;
}

bool OneDMaze::Act(int action)
{
    state = 0;
    reward = 0;

    if (action==0) {
        hidden_state--;
    } else if (action==1) {
        hidden_state++;
    } else {
        Serror("Invalid action %d\n", action);
        exit(-1);
    }

    if (hidden_state>=n_hidden_states) {
        hidden_state = n_hidden_states;
    }
    if (hidden_state==0) {
        hidden_state = 1 + rng->discrete_uniform(n_hidden_states - 1);
        state = 1;
        reward = 1;
    }
    

    
    
}

