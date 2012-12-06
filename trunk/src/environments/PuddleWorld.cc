// -*- Mode: c++ -*-
// copyright (c) 2012 by Nikolaos Tziortziotis
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "PuddleWorld.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include "NormalDistribution.h"

//PuddleWorld::Parameters PuddleWorld::default_parameters = SetDefaultParameters();

PuddleWorld:: PuddleWorld(bool random_parameters)
 : Environment<Vector, int>(2,4)
   //parameters(default_parameters)
{	
    printf("Entering actual constructor\n");

	default_parameters = SetDefaultParameters();
    
	if(random_parameters){
		RandomSourceRNG rng(false);
		parameters.U_POS_X = (0.5 + rng.uniform()) * default_parameters.U_POS_X;
		parameters.L_POS_X = (0.5 + rng.uniform()) * default_parameters.L_POS_X;
		parameters.U_POS_Y = (0.5 + rng.uniform()) * default_parameters.U_POS_Y;
		parameters.L_POS_Y = (0.5 + rng.uniform()) * default_parameters.L_POS_Y;
		parameters.MCNOISE = (0.5 + rng.uniform()) * default_parameters.MCNOISE;

		for(int i=0;i<default_parameters.NUMPUDDLES;i++)
			parameters.RADIUSPUDDLES(i) = (0.5 + rng.uniform()) * default_parameters.RADIUSPUDDLES(i);
		parameters.AGENTSPEED	= (0.5 + rng.uniform()) * default_parameters.AGENTSPEED;
	} else {
        parameters = default_parameters;
    }
	
	state.Resize(n_states);
	state.Clear();
	
	state_upper_bound.Resize(n_states);
	state_lower_bound.Resize(n_states);
	state_upper_bound[0] = parameters.U_POS_X;
	state_upper_bound[1] = parameters.U_POS_Y;
	state_lower_bound[0] = parameters.L_POS_X;
	state_lower_bound[1] = parameters.L_POS_Y;
	
	action_upper_bound.Resize(n_actions);
    action_lower_bound.Resize(n_actions);
	for (uint i=0; i<n_actions; ++i){
		action_lower_bound(i) = 0;
        action_upper_bound(i) = 1;
    }
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
	for(uint i=0; i<n_states; i++){
		state_action_lower_bound(i) = state_lower_bound(i);
		state_action_upper_bound(i) = state_upper_bound(i);
	}
	for(uint i=0; i<n_actions; i++){
		state_action_lower_bound(i + n_states) = action_lower_bound(i);
		state_action_upper_bound(i + n_states) = action_upper_bound(i);
	}
	endsim = false;
}

PuddleWorld::~PuddleWorld(){}

void PuddleWorld::Reset()
{
	state[0] = urandom(parameters.L_POS_Y,parameters.U_POS_Y - 0.1);
	state[1] = urandom(parameters.L_POS_Y,parameters.U_POS_Y - 0.1);
	endsim   = false;
	reward   = -1;
}

bool PuddleWorld::Act(const int& action)
{
	// make sure we tell the guy we have terminated
	if(endsim){
		reward = 0.0;
		return false;
	}
	//run
	Simulate(action);
	
	if(endsim) {
		return false;
	}
	
	return true;
}

void PuddleWorld::Simulate(const int action)
{
	Vector input(2);
	
	switch(action){
    case 0: input[0] = parameters.AGENTSPEED; break;
    case 1: input[0] = -parameters.AGENTSPEED; break;
    case 2: input[1] = parameters.AGENTSPEED; break;
    case 3: input[1] = -parameters.AGENTSPEED; break;
    default: Serror("Undefined action %d\n",action);
	}
	
	state += input;
	
	//We add noise in the transition.
	NormalDistribution R;
	state[0] = state[0] + (R.generate())*parameters.MCNOISE*parameters.AGENTSPEED;
	state[1] = state[1] + (R.generate())*parameters.MCNOISE*parameters.AGENTSPEED;

	if(state[0] > parameters.U_POS_X){
		state[0] = parameters.U_POS_X;
	}
	if(state[0] < parameters.L_POS_X){
		state[0] = parameters.L_POS_X;
	}
	if(state[1] > parameters.U_POS_Y){
		state[1] = parameters.U_POS_Y;
	}
	if(state[1] < parameters.L_POS_Y){
		state[1] = parameters.L_POS_Y;
	}
	
	if(state[0] + state[1] >= 1.9){
		reward = 0.0;
		endsim = true;
	}
	else{
		reward = -1.0;
		for(int i=0;i<parameters.NUMPUDDLES;i++){
			real distance = DistPointToPuddle(i);
			if(distance < parameters.RADIUSPUDDLES[i])
				reward += -400.0*(parameters.RADIUSPUDDLES[i] - distance);
		}
		endsim = false;
	}
	
	return;
}

real PuddleWorld::DistPointToPuddle(const int puddle) const 
{
	Vector P0	= (parameters.U_POS_P).getRow(puddle);
	Vector P1	= parameters.L_POS_P.getRow(puddle);
	
	Vector v	= P1 - P0;
	Vector w	= state - P0;
	
	real c1 = Product(w,v);
	if(c1 <= 0){
		return (state - P0).L2Norm();
	}
	real c2 = Product(v,v);
	if(c2 <= c1){
		return (state - P1).L2Norm();
	}
	
	real c = c1 / c2;
	return (state - (P0 + v*c)).L2Norm();	
}
