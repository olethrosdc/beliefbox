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

#define sqr(x)   ((x)*(x))
Bike::Parameters Bike::default_parameters = 
{
-1.0,
 0.0,
+0.01,
0.0001,
20,
0.01,
(10.0/3.6), 
9.82,
0.3,
0.66,
0.94,
15.0,
1.7,
60.0,
75.0,
0.34,      /* tyre radius */  
1.11,     /* distance between the point where the front and back tyre touch the ground */
3.1415927,
1000,
0,
10
};


Bike::Bike(bool random_parameters)
: Environment<Vector, int>(6, 5), parameters(default_parameters)
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
	Reset();
    state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
	
	sigma_dot = ( ((float) parameters.v) / parameters.R);
	I_bike	= ((13.0/3)*parameters.Mc*parameters.h*parameters.h + parameters.Mp*(parameters.h+parameters.dCM)*(parameters.h+parameters.dCM));
	I_dc	= (parameters.Md*parameters.R*parameters.R);
	I_dv	= ((3.0/2)*parameters.Md*parameters.R*parameters.R); 
	I_dl    = ((1.0/2)*parameters.Md*parameters.R*parameters.R);  
	
	state_lower_bound[0] = -(80.0/180.0)*parameters.pi;
	state_lower_bound[1] = -2.0*parameters.pi;
	state_lower_bound[2] = -parameters.pi/15.0;
	state_lower_bound[3] = -2.0*parameters.pi;
	state_lower_bound[4] = -parameters.pi;
	state_lower_bound[5] = -parameters.pi;

	state_upper_bound[0] = (80.0/180.0)*parameters.pi;
	state_upper_bound[1] = 2.0*parameters.pi;
	state_upper_bound[2] = parameters.pi/15.0;
	state_upper_bound[3] = 2.0*parameters.pi;
	state_upper_bound[4] = parameters.pi;
	state_upper_bound[5] = parameters.pi;
	state_upper_bound.print(stdout);
	state_lower_bound.print(stdout);
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
	omega = omega_dot = omega_d_dot = 0.0;
	theta =  theta_dot = theta_d_dot = 0.0;
	xb = 0.0; yb = parameters.l;
	xf = 0.0; yf = 0.0; // parameters.l;
	psi =  atan((xb-xf)/(yf-yb));
	psi_goal = calc_angle_to_goal_1(xf, xb, yf, yb);
//	psi_goal = parameters.pi / 2.0;
	state[0] = theta;
	state[1] = theta_dot;
	state[2] = omega;
	state[3] = omega_dot;
	state[4] = omega_d_dot;
	state[5] = psi_goal;
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
	static real rCM, rf, rb;
	static real T, d, phi;
	real temp;

	if(action == 0) {
		T = -2.0;
		d = 0.0;
	}
	else if(action == 1) {
		T = 2.0;
		d = 0.0;
	}
	else if(action == 2) {
		T = 0.0; 
		d = -0.02;
	}
	else if(action == 3) {
		T = 0.0;
		d = 0.02;
	}
	else {
		T = 0.0;
		d = 0.0;
	}
	//T = 2.0*(((real)action / 3.0)-1.0);
//	d = 0.02*((real)(action % 3)-1.0);
	d = d + 0.04*(0.5-urandom()); /* Max noise is 2 cm */
//	state.print(stdout);
	if (state[0] == 0.0) {
		rCM = rf = rb = 9999999.0; /* just a large number */
	} else {
		rCM = sqrt(pow(parameters.l-parameters.c,2.0) + parameters.l*parameters.l/(pow(tan(state[0]),2.0)));
		rf = parameters.l / fabs(sin(state[0]));
		rb = parameters.l / fabs(tan(state[0]));
	} /* rCM, rf and rb are always positiv */
	
	/* Main physics eq. in the bicycle model coming here: */
	phi = state[2] + atan(d/parameters.h);
	state[4] = ( parameters.h*parameters.M*parameters.g*sin(phi) 
				   - cos(phi)*(I_dc*sigma_dot*state[1] + sign(state[0])*parameters.v*parameters.v*
							   (parameters.Md*parameters.R*(1.0/rf + 1.0/rb) +  parameters.M*parameters.h/rCM))
				   ) / I_bike;
	theta_d_dot =  (T - I_dv*state[3]*sigma_dot) /  I_dl;
	
	/*--- Eulers method ---*/
	state[3] += state[4] * parameters.dt;
	state[2] += state[3] * parameters.dt;
	state[1] += theta_d_dot * parameters.dt;
	state[0] += state[1] * parameters.dt;
	
	if (fabs(state[0]) > 1.3963) { /* handlebars cannot turn more than 
	 80 degrees */
		state[0] = sign(state[0]) * 1.3963;
	}
	
	/* New position of front tyre */
	temp = parameters.v*parameters.dt/(2.0*rf);                             
	if (temp > 1) temp = sign(psi + state[0]) * parameters.pi/2.0;
	else temp = sign(psi + state[0]) * asin(temp); 
	xf += parameters.v * parameters.dt * (-sin(psi + state[0] + temp));
	yf += parameters.v * parameters.dt * cos(psi + state[0] + temp);
	
	/* New position of back tyre */
	temp = parameters.v*parameters.dt/(2.0*rb);               
	if (temp > 1) temp = sign(psi) * parameters.pi/2.0;
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
	if ((xf == xb) && (temp < 0.0)) psi = parameters.pi;
	else {
		if (temp > 0.0) psi = atan((xb-xf)/temp);
		else psi = sign(xb-xf)*(parameters.pi/2.0) - atan(temp/(xb-xf));
	}
	
	state[5] = calc_angle_to_goal(xf, xb, yf, yb);
	/*-- Calculation of the reinforcement  signal --*/
	if (fabs(state[2]) > (parameters.pi/15.0)) { /* the bike has fallen over */
		reward = parameters.R1;
		endsim = true;
		/* a good place to print some info to a file or the screen */  
	} else { 
		temp = calc_dist_to_goal(xf, xb, yf, yb);
		
		if (temp < 1e-3) reward = parameters.R3;
		else reward = (0.95 - sqr(state[5])) * parameters.R_FACTOR; 
		endsim = false;
	}
	state[5] = acos(calc_angle_to_goal_1(xf, xb, yf, yb));

//	printf("reward = %f\n",reward);
}

real Bike::calc_dist_to_goal(real xf, real xb, real yf, real yb)
{
	real temp;
	
	temp = sqrt(std::max(0.0, (parameters.x_goal-xf)*(parameters.x_goal-xf) + (parameters.y_goal-yf)*(parameters.y_goal-yf) 
					- parameters.radius_goal*parameters.radius_goal)); 
	return(temp);
}


real Bike::calc_angle_to_goal(real xf, real xb, real yf, real yb)
{
	real temp, scalar, tvaer;
	
	temp = (xf-xb)*(parameters.x_goal-xf) + (yf-yb)*(parameters.y_goal-yf); 
//	printf("temp = %f\n",temp);
	scalar =  temp / (parameters.l * sqrt(sqr(parameters.x_goal-xf)+sqr(parameters.y_goal-yf)));
//	printf("scalar = %f\n",scalar);
	tvaer = (-yf+yb)*(parameters.x_goal-xf) + (xf-xb)*(parameters.y_goal-yf); 
//	printf("tvaer = %f\n",tvaer);
	if (tvaer <= 0.0) temp = scalar - 1.0;
	else temp = fabs(scalar - 1.0);
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

real Bike::calc_angle_to_goal_1(real xf, real xb, real yf, real yb)
{
	real temp	= (xf - xb)*(parameters.x_goal - xf) + (yf-yb)*(parameters.y_goal - yf);
	real scalar = temp / (parameters.l * sqrt(sqr(parameters.x_goal-xf)+sqr(parameters.y_goal-yf)));
	return(scalar);
}


real Bike::sign(const real& num) {
	if(num == 0) {
		return 0.0;
	} else if(num > 0.0) {
		return 1.0;
	} else {
		return -1.0;
	}
	
}

