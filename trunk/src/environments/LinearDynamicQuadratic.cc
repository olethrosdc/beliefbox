/*
 *  LinearDynamicQuadratic.cc
 *  
 *
 *  Created by Poseidon on 20/12/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "LinearDynamicQuadratic.h"
#include "Random.h"
#include "RandomSourceRNG.h"

LinearDynamicQuadratic::Parameters LinearDynamicQuadratic::default_parameters = 
{
	4.0,		//Upper bound on position
	-4.0,		//Lower bound on position
	0.1,		//Upper bound on velocity
	-0.1,		//Lower bound on velocity
	0.1,		//Contribution on input
};

LinearDynamicQuadratic::LinearDynamicQuadratic(bool random_parameters)
: Environment<Vector, int>(2,3), parameters(default_parameters)
{
	if(random_parameters) {
		RandomSourceRNG rng(false);
		parameters.U_POS = (0.5 + rng.uniform()) * default_parameters.U_POS;
        parameters.L_POS = (0.5 + rng.uniform()) * default_parameters.L_POS;
        parameters.U_VEL = (0.5 + rng.uniform()) * default_parameters.U_VEL;
        parameters.L_VEL = (0.5 + rng.uniform()) * default_parameters.L_VEL;
        parameters.INPUT = (0.5 + rng.uniform()) * default_parameters.INPUT;
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

LinearDynamicQuadratic::~LinearDynamicQuadratic()
{
	// nothing to do
}

void LinearDynamicQuadratic::Reset()
{
	state[0] = urandom(parameters.L_POS, parameters.U_POS);
    state[1] = urandom(parameters.L_VEL, parameters.U_VEL);
	// state[0] = -0.5;
	//	state[1] = 0.0;
	
	endsim = false;
    reward = 0.0;
}

bool LinearDynamicQuadratic::Act(const int& action)
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

void LinearDynamicQuadratic::Simulate(const int& action)
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
    
	real vel = state[1];
    state[1] = state[1] + parameters.INPUT*input;
    if (state[1] > parameters.U_VEL) {
        state[1] = parameters.U_VEL;
    }
    if (state[1] < parameters.L_VEL) {
        state[1] = parameters.L_VEL;
    }
	
	vel = ((vel + state[1])*parameters.INPUT)/2;
	
    state[0] = state[0] + vel;
    if (state[0] > parameters.U_POS) {
        state[0] = parameters.U_POS - 0.01;
		state[1] = 0.0;
    }
    if (state[0] < parameters.L_POS) {
        state[0] = parameters.L_POS + 0.01;
        state[1] = 0.0;
    }
	
    if (state[0] > -0.5 && state[0] < 0.5 ) {
        reward = 0.0;
        endsim = true;
    } else {
        reward = -1;
        endsim = false;
    }
	
    return;
	
}
