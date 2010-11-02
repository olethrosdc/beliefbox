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

#include "ContinuousChain.h"
#include "Random.h"

ContinuousChain::ContinuousChain() : Environment<Vector, int> (1, 3)
{
    state.Resize(1);
    endsim = false;
}

ContinuousChain::~ContinuousChain()
{
    // nothing to do
}
void ContinuousChain::Reset()
{
    state[0] = 0.0;// urandom(-, U_POS);
    endsim = false;
    reward = 0.0;
}

bool ContinuousChain::Act(int action)
{
    // make sure we tell the guy we have terminated
    if (endsim) {
        reward = 0.0;
        return false;
    }
    
    // run
    Simulate(action);
    return !endsim;
}

void ContinuousChain::Simulate(int action)
{
    real input=0.0;

    switch (action){
    case 0: input = -1.0; break;
    case 1: input = 0.0; break;
    case 2: input = 1.0; break;
    default: Serror("Undefined action %d\n", action);
    }

    
    state[0] += urandom()*.5*input;

    reward = -0.1;
    if (state[0] >= 1.0) {
        reward = 1.0;
        state[0] = 1.0;
    }

    if (state[0] <= -1.0) {
        reward = -1.0;
        state[0] = -1.0;
    }



    return;
  
}

