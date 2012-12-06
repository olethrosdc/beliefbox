// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// Adapted from code by Michail G. Lagoudakis, copyright (c) 2003-2008 
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
#include "RandomSourceRNG.h"

MountainCar::Parameters MountainCar::default_parameters = 
    {
        0.5,            // Upper bound on position
        -1.2,           // Lower bound on position
        0.07,           // Upper bound on velocity
        -0.07,          // Lower bound on velocity
        0.001,          // contribution of input
        0.0025,         // contribution of gravity
        0.2             // input noise
    };
    
    MountainCar::MountainCar(bool random_parameters)
        : Environment<Vector, int>(2, 3),
          parameters(default_parameters)
                                 
{
    if (random_parameters) {
        RandomSourceRNG rng(false);
        parameters.U_POS = (0.5 + rng.uniform()) * default_parameters.U_POS;
        parameters.L_POS = (0.5 + rng.uniform()) * default_parameters.L_POS;
        parameters.U_VEL = (0.5 + rng.uniform()) * default_parameters.U_VEL;
        parameters.L_VEL = (0.5 + rng.uniform()) * default_parameters.L_VEL;
        parameters.INPUT = (0.5 + rng.uniform()) * default_parameters.INPUT;
        parameters.GRAVITY = (0.5 + rng.uniform()) * default_parameters.GRAVITY;
        parameters.MCNOISE = (0.5 + rng.uniform()) * default_parameters.MCNOISE;
    }
    state.Resize(n_states);
	state.Clear();

    state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
    state_upper_bound[0] = parameters.U_POS;
    state_upper_bound[1] = parameters.U_VEL;
    state_lower_bound[0] = parameters.L_POS;
    state_lower_bound[1] = parameters.L_VEL;
   
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

MountainCar::~MountainCar()
{
    // nothing to do
}


void MountainCar::Reset()
{
    state[0] = -0.5;// urandom(-, U_POS);
    state[1] = 0.0;//urandom(L_VEL, U_VEL);
    endsim = false;
    reward = 0.0;
}
bool MountainCar::Act(const int& action)
{
    // make sure we tell the guy we have terminated
    if (endsim) {
        reward = 0.0;
        return false;
    }
    
    // run
    Simulate(action);

	if (endsim) {
		return false;
	}
    return true;
}

void MountainCar::Simulate(const int action)
{
    real input=0.0;

    switch (action){
    case 0: input = -1.0; break;
    case 1: input = 0.0; break;
    case 2: input = 1.0; break;
    default: Serror("Undefined action %d\n", action);
    }

    real noise = urandom(-parameters.MCNOISE, parameters.MCNOISE);
    input += noise;
    
    state[1] = state[1] + parameters.INPUT*input - parameters.GRAVITY*cos(3.0*state[0]);
    if (state[1] > parameters.U_VEL) {
        state[1] = parameters.U_VEL;
    }
    if (state[1] < parameters.L_VEL) {
        state[1] = parameters.L_VEL;
    }

    state[0] = state[0] + state[1];
    if (state[0] > parameters.U_POS) {
        state[0] = parameters.U_POS;
    }
    if (state[0] < parameters.L_POS) {
        state[0] = parameters.L_POS + 0.01;
        state[1] = 0.01;
    }
    //printf ("S: %f %f\n", state[0], state[1]);
    if (state[0] == parameters.U_POS) {
        reward = 0.0;
        endsim = true;
    } else {
        reward = -1;
        endsim = false;
    }
  

    return;
  
}

