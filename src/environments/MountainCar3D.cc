// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MountainCar3D.h"
#include "Random.h"

real MountainCar3D::U_POS = 0.6;                // Upper bound on position
real MountainCar3D::L_POS = -1.2;             // Lower bound on position
real MountainCar3D::U_VEL = 0.07;               // Upper bound on velocity
real MountainCar3D::L_VEL = -0.07;            // Lower bound on velocity
real MountainCar3D::INPUT = 0.001;           // contribution of input
real MountainCar3D::GRAVITY = 0.0025;          // contribution of gravity
real MountainCar3D::MCNOISE = 0.2;

MountainCar3D::MountainCar3D() : Environment<Vector, int>(4, 4)
{
    state.Resize(n_states);
   
	state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
    state_upper_bound[0] = U_POS;
    state_upper_bound[1] = U_POS;
    state_upper_bound[2] = U_VEL;
    state_upper_bound[3] = U_VEL;
    state_lower_bound[0] = L_POS;
    state_lower_bound[1] = L_POS;
    state_lower_bound[2] = L_VEL;
    state_lower_bound[3] = L_VEL;
   
	action_upper_bound.Resize(n_actions);
	action_lower_bound.Resize(n_actions);
	for (uint i=0; i<n_actions; ++i) {
		action_lower_bound(i) = 0;
		action_upper_bound(i) = 1;
	}
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
	for (uint i=0; i<n_states; ++i) {
		state_action_lower_bound(i) = state_lower_bound(i);
		state_action_upper_bound(i) = state_upper_bound(i);
	}
	for (uint i=0; i<n_actions; ++i) {
		state_action_lower_bound(i + n_states) = action_lower_bound(i);
		state_action_upper_bound(i + n_states) = action_upper_bound(i);
	}
    endsim = false;
}

MountainCar3D::~MountainCar3D()
{
    // nothing to do
}
void MountainCar3D::Reset()
{
    state[0] = 0.0;// urandom(-, U_POS);
    state[1] = 0.0;// urandom(-, U_POS); 
    state[2] = 0.0;//urandom(L_VEL, U_VEL);
    state[3] = 0.0;//urandom(L_VEL, U_VEL);
    endsim = false;
    reward = -1;
}
bool MountainCar3D::Act(int action)
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

void MountainCar3D::Simulate(int action)
{
    real input_x = 0;
    real input_y = 0;
    switch (action){
    case 0: input_x = 1.0; break;
    case 1: input_x = -1.0; break;
    case 2: input_y = 1.0; break;
    case 3: input_y = -1.0; break;
    default: Serror("Undefined action %d\n", action);
    }
    
    input_x += urandom(-MCNOISE, MCNOISE);
    input_y += urandom(-MCNOISE, MCNOISE);
    
    state[2] += INPUT*input_x - GRAVITY*cos(3.0*state[0]);
    state[3] += INPUT*input_y - GRAVITY*cos(3.0*state[1]);

    state[0] += state[2];
    state[1] += state[3];
        
    for (uint i=0; i<n_states; ++i) {
        if (state[i] > state_upper_bound[i]) {
            state[i] = state_upper_bound[i];
        }
        if (state[i] < state_lower_bound[i]) {
            state[i] = state_lower_bound[i];
        }
    }
    
    reward = -1;
    endsim = false;
    if ((state[0] >= state_upper_bound[0] - 10e-6)
        && (state[1] >= state_upper_bound[1] - 10e-6)) {
        reward = 1;
        endsim = true;
    }
    return;
  
}

