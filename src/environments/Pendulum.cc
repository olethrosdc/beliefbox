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

#include "Pendulum.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include "MersenneTwister.h"

Pendulum::Parameters Pendulum::default_parameters = 
    { 
        2.0, // pendulum mass
        8.0, // cart mass
        0.5, // pendulum length
        9.8, // gravity
        0.01, // max noise
        0.01 // interval
    };

Pendulum::Pendulum(bool random_parameters)
    : Environment<Vector, int>(2, 3),
      parameters(default_parameters),
      CCa (1.0 / (parameters.pendulum_mass + parameters.cart_mass))
{
    if (random_parameters) {
        //RandomSourceRNG rng(false);
		//MersenneTwisterRNG rng;
		//rng.manualSeed(12315);
        parameters.pendulum_mass = (0.5 + urandom()) * default_parameters.pendulum_mass;
        parameters.cart_mass = (0.5 + urandom()) * default_parameters.cart_mass;
        parameters.pendulum_length = (0.5 + urandom()) * default_parameters.pendulum_length;
        parameters.gravity = (0.5 + urandom()) * 
default_parameters.gravity;
        parameters.max_noise = (0.5 + urandom()) * default_parameters.max_noise;
        parameters.Dt = (0.5 + urandom()) * default_parameters.Dt;
    }

    state.Resize(2);
    state.Clear();
    state_upper_bound.Resize(2);
    state_lower_bound.Resize(2);
    state_upper_bound[0] = 1.6;
    state_upper_bound[1] = -1.5;//10;
    state_lower_bound[0] = -1.6;
    state_lower_bound[1] = 1.5;//-10;

	action_upper_bound.Resize(n_actions);
	action_lower_bound.Resize(n_actions);
	action_upper_bound += 1;
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
    
    reward = 0;

    endsim = false;
}

Pendulum::~Pendulum()
{
    // nothing to do
}

void Pendulum::Reset()
{
    reward = 1.0;
#if 1
    // Theta
    state[0] =  urandom(-0.01, 0.01);
    // dTheta/dt
    state[1] = urandom(-0.001, 0.001);
#else
	for (int i=0; i<2; ++i) {
		state[i] = urandom(state_lower_bound[i], state_upper_bound[i]);
	}
#endif
    endsim = false;
}

void Pendulum::penddot(Vector& xdot, real u, Vector& x)
{
    // Nonlinear model 
     
    double cx = cos(x[0]);
    real dtheta2 =x[1]*x[1];
    xdot[0] = x[1]; 
    xdot[1] = (parameters.gravity * sin(x[0]) - 
               0.5*CCa * parameters.pendulum_mass * parameters.pendulum_length * dtheta2 * sin(2.0*x[0]) -  CCa * cos(x[0]) * u ) / 
        ( 4.0/3.0*parameters.pendulum_length - CCa*parameters.pendulum_mass*parameters.pendulum_length*cx*cx ); 
   
}

bool Pendulum::Act(const int& action)
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

void Pendulum::Simulate(const int action)
{
    Vector xdot(2);
    real input=0.0, noise, t;

    //printf ("# s: %f %f, a: %d\n", state[0], state[1], action);
    switch(action) {
    case 0: input = -50.0; break;
    case 1: input = 0.0; break;
    case 2:  input = +50.0; break;
    }

    noise = urandom(-parameters.max_noise, parameters.max_noise);
    input += noise;

    // Simulate for 0.1 seconds
    for (t=0.0; t<=0.1; t+=parameters.Dt) {

        penddot(xdot, input, state);
    
        state[0] += xdot[0] * parameters.Dt;
        state[1] += xdot[1] * parameters.Dt;

    }
  
    if (fabs(state[0]) > M_PI/2.0) {
        reward = -1.0;
        endsim = true;
    } else {
        reward = 0.0;
        endsim = false;
    }
}
