// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// copyright (c) 2003-2008 Michail G. Lagoudakis
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MountainCar.h"
#include "Random.h"

real MountainCar::U_POS = 0.5;                // Upper bound on position
real MountainCar::L_POS = -1.2;             // Lower bound on position
real MountainCar::U_VEL = 0.07;               // Upper bound on velocity
real MountainCar::L_VEL = -0.07;            // Lower bound on velocity
real MountainCar::INPUT = 0.001;           // contribution of input
real MountainCar::GRAVITY = 0.0025;          // contribution of gravity
real MountainCar::MCNOISE = 0.2;

MountainCar::MountainCar() : Environment<Vector, int>(2, 3)
{
    state.Resize(n_states);
    endsim = false;
}

MountainCar::~MountainCar()
{
    // nothing to do
}
void MountainCar::Reset()
{
    state[0] = 0.0;// urandom(-, U_POS);
    state[1] = 0.0;//urandom(L_VEL, U_VEL);
    endsim = false;
    reward = 0.0;
}
bool MountainCar::Act(int action)
{
    // make sure we tell the guy we have terminated
    if (endsim) {
        reward = 0.0;
        return false;
    }
    
    // run
    Simulate(action);
    return true;
}

void MountainCar::Simulate(int action)
{
    real input=0.0;

    switch (action){
    case 0: input = -1.0; break;
    case 1: input = 0.0; break;
    case 2: input = 1.0; break;
    default: Serror("Undefined action %d\n", action);
    }

    real noise = urandom(-MCNOISE, MCNOISE);
    input += noise;
    
    state[1] = state[1] + INPUT*input - GRAVITY*cos(3.0*state[0]);
    if (state[1] > U_VEL) {
        state[1] = U_VEL;
    }
    if (state[1] < L_VEL) {
        state[1] = L_VEL;
    }

    state[0] = state[0] + state[1];
    if (state[0] > U_POS) {
        state[0] = U_POS;
    }
    if (state[0] < L_POS) {
        state[0] = L_POS + 0.01;
        state[1] = 0.01;
    }
        //printf ("S: %f %f\n", state[0], state[1]);
    if (state[0] == U_POS) {
        reward = 0.0;
        endsim = true;
    } else {
        reward = -0.01;
        endsim = false;
    }
  

    return;
  
}

