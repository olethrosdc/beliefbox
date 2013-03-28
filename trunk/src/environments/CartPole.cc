// -*- Mode: c++ -*-
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@cs.uoi.gr>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "CartPole.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include "MersenneTwister.h"

CartPole::Parameters CartPole::default_parameters = 
{ 
9.8, 
1.0, 
0.1, 
0.5, 
10.0,
0.02,
0.0
};

const real CartPole::FOURTHIRDS = 4.0 / 3.0;



CartPole::CartPole(bool random_parameters)
: Environment<Vector, int>(4, 3),
parameters(default_parameters)
{
    if (random_parameters) {
        //RandomSourceRNG rng(false);
		//MersenneTwisterRNG rng;
		//rng.manualSeed(12315);
        parameters.GRAVITY = (0.5 + urandom()) * default_parameters.GRAVITY;
        parameters.MASSCART = (0.5 + urandom()) * default_parameters.MASSCART;
        parameters.MASSPOLE = (0.5 + urandom()) * default_parameters.MASSPOLE;
        parameters.LENGTH = (0.5 + urandom()) * default_parameters.LENGTH;
//        parameters.max_noise = (0.5 + urandom()) * default_parameters.max_noise;
		parameters.FORCE_MAG = (0.5 + urandom()) * default_parameters.FORCE_MAG;
        parameters.TAU = (0.5 + urandom()) * default_parameters.TAU;
		parameters.noise = (0.5 + urandom()) * default_parameters.noise;
    }
	TOTAL_MASS = (parameters.MASSPOLE + parameters.MASSCART);
	POLEMASS_LENGTH = (parameters.MASSPOLE * parameters.LENGTH);
	
	
	state.Resize(n_states);
    state.Clear();
    state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
	state_lower_bound[0] = -2.4;
    state_upper_bound[0] = 2.4;
    state_lower_bound[1] = -6.0;//-10;
    state_upper_bound[1] = 6.0;//10;
    state_lower_bound[2] = -0.2094;
    state_upper_bound[2] = 0.2094;
    state_lower_bound[3] = -6.0;
    state_upper_bound[3] = 6.0; 
	
	action_upper_bound.Resize(n_actions);
	action_lower_bound.Resize(n_actions);
	action_upper_bound += 1;
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
    
    reward = 0;
	
    endsim = false;
}

CartPole::~CartPole()
{
    // nothing to do
}

void CartPole::Reset()
{
    reward = 1.0;
#if 1
    /// Cart position
    state[0] = 0.0; // urandom(-0.01, 0.01);
	/// Cart velocity
    state[1] = 0.0;
	// Theta
    state[2] = 0.0;
	// dTheta/dt
    state[3] = urandom(-0.001, 0.001);
#else
	for (int i=0; i<4; ++i) {
		state[i] = urandom(state_lower_bound[i], state_upper_bound[i]);
	}
#endif
    endsim = false;
}

bool CartPole::Act(const int& action)
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

void CartPole::Simulate(const int action)
{
    real xacc;
	real thetaacc;
	real costheta;
	real sintheta;
	real temp;
	real force = 0.0;;
	
	switch(action) {
		case 0: force = -parameters.FORCE_MAG; break;
		case 1: force = 0.0; break;
		case 2: force = parameters.FORCE_MAG; break;
    }
	
	//Noise of 1.0 means possibly full opposite action
	real thisNoise=2.0*parameters.noise*parameters.FORCE_MAG*(urandom()-0.5);
	
	force+=thisNoise;
	
	costheta = cos(state[2]);
	sintheta = sin(state[2]);
	
	temp = (force + POLEMASS_LENGTH * state[3]* state[3]* sintheta) / TOTAL_MASS;
	thetaacc = (parameters.GRAVITY * sintheta - costheta * temp) / (parameters.LENGTH * (FOURTHIRDS - parameters.MASSPOLE * costheta * costheta / TOTAL_MASS));

	xacc = temp - POLEMASS_LENGTH * thetaacc * costheta / TOTAL_MASS;
	/*** Update the four state variables, using Euler's method. ***/
	state[0] += parameters.TAU * state[1];
	state[1] += parameters.TAU * xacc;
	state[2] += parameters.TAU * state[3];
	state[3] += parameters.TAU * thetaacc;
	
	/**These probably never happen because the pole would crash **/
	while (state[2] >= M_PI) {
		state[2] -= 2.0 * M_PI;
	}
	while (state[2] < -M_PI) {
		state[2] += 2.0 * M_PI;
	}
	
	if (state[0] < state_lower_bound[0] || state[0] > state_upper_bound[0] || state[2] < state_lower_bound[2] || state[2] > state_upper_bound[2]) {
		endsim = true;
		reward = -1.0;
	} 
	else {
		endsim = false;
		reward = 1.0;
	}
	
}

