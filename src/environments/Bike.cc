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

#include "Bike.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include "MersenneTwister.h"

Bike::Parameters Bike::default_parameters = 
	{
		-1.0, // r1
		0.0, // r2
		+0.01, // r3
		0.0001, // r_factor
		20, // n0
		0.01, // dt
		(10.0/3.6), // v 
		9.82, // g
		0.3, // dC
		0.66, // c
		0.94, // h
		15.0, // Mc
		1.7, // Md
		60.0, // Mp
		75.0, // M
		0.34,      /* tyre radius */  
		1.11,     /* distance between the point where the front and back tyre touch the ground */
		3.1415927, // pi
		1000, // x_goal
		0, // y_goal
		10, // r_goal
		0.04 // max_noise
	};


Bike::Bike(bool random_parameters)
  : Environment<Vector, int>(5, 9), parameters(default_parameters)
{
    if (random_parameters) {
        //RandomSourceRNG rng(false);
		//MersenneTwisterRNG rng;
		//rng.manualSeed(12315);
        parameters.R1 = (0.5 + urandom()) * default_parameters.R1;
        parameters.R2 = (0.5 + urandom()) * default_parameters.R2;
        parameters.R3 = (0.5 + urandom()) * default_parameters.R3;
        parameters.R_FACTOR = (0.5 + urandom()) * default_parameters.R_FACTOR;
        parameters.N0_STATES2 = (0.5 + urandom()) * default_parameters.N0_STATES2;
		parameters.dt = (0.5 + urandom()) * default_parameters.dt;
		parameters.v = (0.5 + urandom()) * default_parameters.v;
        parameters.g = (0.5 + urandom()) * default_parameters.g;
        parameters.dCM = (0.5 + urandom()) * default_parameters.dCM;
		parameters.c = (0.5 + urandom()) * default_parameters.c;
		parameters.h = (0.5 + urandom()) * default_parameters.h;
        parameters.Mc = (0.5 + urandom()) * default_parameters.Mc;
        parameters.Md = (0.5 + urandom()) * default_parameters.Md;
        parameters.Mp = (0.5 + urandom()) * default_parameters.Mp;
		parameters.M = (0.5 + urandom()) * default_parameters.M;
        parameters.R = (0.5 + urandom()) * default_parameters.R;
        parameters.l = (0.5 + urandom()) * default_parameters.l;
		parameters.x_goal = (0.5 + urandom()) * default_parameters.x_goal;
        parameters.y_goal = (0.5 + urandom()) * default_parameters.y_goal;
        parameters.radius_goal = (0.5 + urandom()) * default_parameters.radius_goal;
    }
	
    state.Resize(n_states);
    state.Clear();
    state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
	
	sigma_dot = ( ((float) parameters.v) / parameters.R);
	I_bike	= ((13.0/3)*parameters.Mc*parameters.h*parameters.h + parameters.Mp*(parameters.h+parameters.dCM)*(parameters.h+parameters.dCM));
	I_dc	= (parameters.Md*parameters.R*parameters.R);
	I_dv	= ((3.0/2)*parameters.Md*parameters.R*parameters.R); 
	I_dl    = ((1.0/2)*parameters.Md*parameters.R*parameters.R);  
	
	state_lower_bound[0] = -1.5; 
	state_lower_bound[1] = -3;
	state_lower_bound[2] = -0.2;
	state_lower_bound[3] = -0.5;
	state_lower_bound[4] = -2.0;
	//	state_lower_bound[5] = ;

	state_upper_bound[0] = 1.5;
	state_upper_bound[1] = 3;
	state_upper_bound[2] = 0.2;
	state_upper_bound[3] = 0.5;
	state_upper_bound[4] = 2.0;
	//	state_upper_bound[5] = 

	action_upper_bound.Resize(n_actions);
	action_lower_bound.Resize(n_actions);
	action_upper_bound += 1;
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
    
    reward = 0;
	
    endsim = false;
}

Bike::~Bike()
{
    // nothing to do
}

void Bike::Reset()
{
	reward = -1.0;
	omega = omega_dot = omega_d_dot = 0;
	theta =  theta_dot = theta_d_dot = 0;
	xb = 0; yb = 0;
	xf = 0; yf = parameters.l;
	psi =  atan((xb-xf)/(yf-yb));
	psi_goal = calc_angle_to_goal(xf, xb, yf, yb);
	state[0] = theta;
	state[1] = theta_dot;
	state[2] = omega;
	state[3] = omega_dot;
	state[4] = omega_d_dot;
	//	state[5] = psi_goal;
    endsim = false;
}

bool Bike::Act(const int& action)
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

void Bike::Simulate(const int action)
{
	static double rCM, rf, rb;
	static float T, d, phi;
	real temp;
	
	T = 2*((action / 3)-1);
	d = 0.02*((action % 3)-1);
	d = d + parameters.max_noise*(0.5-random()); /* Max noise is 2 cm */
	
	if (theta == 0) {
		rCM = rf = rb = 9999999; /* just a large number */
	} else {
		rCM = sqrt(pow(parameters.l-parameters.c,2) + parameters.l*parameters.l/(pow(tan(theta),2)));
		rf = parameters.l / fabs(sin(theta));
		rb = parameters.l / fabs(tan(theta));
	} /* rCM, rf and rb are always positiv */
	
	/* Main physics eq. in the bicycle model coming here: */
	phi = omega + atan(d/parameters.h);
	omega_d_dot = ( parameters.h*parameters.M*parameters.g*sin(phi) 
					- cos(phi)*(I_dc*sigma_dot*theta_dot + sign(theta)*parameters.v*parameters.v*
								(parameters.Md*parameters.R*(1.0/rf + 1.0/rb) +  parameters.M*parameters.h/rCM))
					) / I_bike;
	theta_d_dot =  (T - I_dv*omega_dot*sigma_dot) /  I_dl;
	
	/*--- Eulers method ---*/
	omega_dot += omega_d_dot * parameters.dt;
	omega += omega_dot * parameters.dt;
	theta_dot += theta_d_dot * parameters.dt;
	theta += theta_dot * parameters.dt;
	
	if (fabs(theta) > 1.3963) { /* handlebars cannot turn more than 
								   80 degrees */
		theta = sign(theta) * 1.3963;
	}
	
	/* New position of front tyre */
	temp = parameters.v*parameters.dt/(2*rf);                             
	if (temp > 1) temp = sign(psi + theta) * parameters.pi/2;
	else temp = sign(psi + theta) * asin(temp); 
	xf += parameters.v * parameters.dt * (-sin(psi + theta + temp));
	yf += parameters.v * parameters.dt * cos(psi + theta + temp);
	
	/* New position of back tyre */
	temp = parameters.v*parameters.dt/(2*rb);               
	if (temp > 1) temp = sign(psi) * parameters.pi/2;
	else temp = sign(psi) * asin(temp); 
	xb += parameters.v * parameters.dt * (-sin(psi + temp));
	yb += parameters.v * parameters.dt * (cos(psi + temp));
	
	/* Round off errors accumulate so the length of the bike changes over many
	   iterations. The following take care of that: */
	temp = sqrt((xf-xb)*(xf-xb)+(yf-yb)*(yf-yb));
	if (fabs(temp - parameters.l) > 0.01) {
		xb += (xb-xf)*(parameters.l-temp)/temp;
		yb += (yb-yf)*(parameters.l-temp)/temp;
	}
	
	temp = yf - yb;
	if ((xf == xb) && (temp < 0)) psi = parameters.pi;
	else {
		if (temp > 0) psi = atan((xb-xf)/temp);
		else psi = sign(xb-xf)*(parameters.pi/2) - atan(temp/(xb-xf));
	}
	
	psi_goal = calc_angle_to_goal(xf, xb, yf, yb);
	
	/*-- Calculation of the reinforcement  signal --*/
	if (fabs(omega) > (parameters.pi/15)) { /* the bike has fallen over */
		reward = parameters.R1;
		endsim = true;
		/* a good place to print some info to a file or the screen */  
	} else { 
		temp = calc_dist_to_goal(xf, xb, yf, yb);
		
		if (temp < 1e-3) reward = parameters.R3;
		else reward = (0.95 - sqrt(psi_goal)) * parameters.R_FACTOR; 
		endsim = false;
	}
	
	state[0] = theta;
	state[1] = theta_dot;
	state[2] = omega;
	state[3] = omega_dot;
	state[4] = omega_d_dot;
	//	state[5] = psi_goal;
}

float Bike::calc_dist_to_goal(float xf, float xb, float yf, float yb)
{
	float temp;
	
	temp = sqrt(std::max(0.0, (parameters.x_goal-xf)*(parameters.x_goal-xf) + (parameters.y_goal-yf)*(parameters.y_goal-yf) 
						 - parameters.radius_goal*parameters.radius_goal)); 
	return(temp);
}


float Bike::calc_angle_to_goal(float xf, float xb, float yf, float yb)
{
	real temp, scalar, tvaer;
	
	temp = (xf-xb)*(parameters.x_goal-xf) + (yf-yb)*(parameters.y_goal-yf); 
	scalar =  temp / (parameters.l * sqrt(sqrt(parameters.x_goal-xf)+sqrt(parameters.y_goal-yf)));
	tvaer = (-yf+yb)*(parameters.x_goal-xf) + (xf-xb)*(parameters.y_goal-yf); 
	
	if (tvaer <= 0) temp = scalar - 1;
	else temp = fabs(scalar - 1);
	
	/* These angles are neither in degrees nor radians, but something
	   strange invented in order to save CPU-time. The measure is arranged the
	   same way as radians, but with a slightly different negative factor. 
	 
	   Say, the goal is to the east.
	   If the agent rides to the east then  temp = 0
	   - " -          - " -   north              = -1
	   - " -                  west               = -2 or 2
	   - " -                  south              =  1 */ 
	
	return(temp);
}
