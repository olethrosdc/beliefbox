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

#include "DiscretisedContinuousChain.h"
#include "Random.h"

DiscretisedContinuousChain::DiscretisedContinuousChain(int n_states_) : Environment<int, int> (n_states_, 3)
{
    endsim = false;
}

DiscretisedContinuousChain::~DiscretisedContinuousChain()
{
    // nothing to do
}
void DiscretisedContinuousChain::Reset()
{
    position = 0.0;// urandom(-, U_POS);
    endsim = false;
    reward = 0.0;
    DiscretiseState();
}

bool DiscretisedContinuousChain::Act(int action)
{
    // make sure we tell the guy we have terminated
    if (endsim) {
        reward = 0.0;
        state = 0;
        return false;
    }
    
    // run
    Simulate(action);
    DiscretiseState();
    return !endsim;
}

void DiscretisedContinuousChain::Simulate(int action)
{
    real input=0.0;

    switch (action){
    case 0: input = -1.0; break;
    case 1: input = 0.0; break;
    case 2: input = 1.0; break;
    default: Serror("Undefined action %d\n", action);
    }

    
    position += 0.1*input;

    reward = -0.1;
    if (position >= 1.0) {
        reward = 1.0;
        position = 1.0;
        endsim = true;
    }

    if (position <= -1.0) {
        reward = -1.0;
    }

    return;
  
}

