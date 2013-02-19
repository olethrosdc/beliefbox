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
 

#include "Acrobot.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include "MersenneTwister.h"
  
Acrobot::Parameters Acrobot::default_parameters = 
{
M_PI,
M_PI,
4*M_PI,
9*M_PI,
1.0,
1.0,
1.0,
1.0,
0.5,
0.5,
1.0,
1.0,
9.8,
0.05,
1.0,
0.0
};

Acrobot::Acrobot(bool random_parameters)
: Environment<Vector, int>(4, 3), parameters(default_parameters)
{
    if (random_parameters) {
        //RandomSourceRNG rng(false);
		//MersenneTwisterRNG rng;
		//rng.manualSeed(12315);
        parameters.m1 = (0.5 + urandom()) * default_parameters.m1;
        parameters.m2 = (0.5 + urandom()) * default_parameters.m2;
        parameters.l1 = (0.5 + urandom()) * default_parameters.l1;
        parameters.l2 = (0.5 + urandom()) * default_parameters.l2;
        parameters.lc1 = (0.5 + urandom()) * default_parameters.lc1;
		parameters.lc2 = (0.5 + urandom()) * default_parameters.lc2;
		parameters.I1 = (0.5 + urandom()) * default_parameters.I1;
        parameters.I2 = (0.5 + urandom()) * default_parameters.I2;
        parameters.g = (0.5 + urandom()) * default_parameters.g;
		parameters.dt = (0.5 + urandom()) * default_parameters.dt;
		parameters.acrobotGoalPosition = (0.5 + urandom()) * default_parameters.acrobotGoalPosition;
	}
	
    state.Resize(n_states);
    state.Clear();
	Reset();
    state_upper_bound.Resize(n_states);
    state_lower_bound.Resize(n_states);
	
	state_upper_bound[0] = parameters.maxTheta1;
	state_upper_bound[1] = parameters.maxTheta2;
	state_upper_bound[2] = parameters.maxTheta1Dot;
	state_upper_bound[3] = parameters.maxTheta2Dot;
	action_upper_bound.Resize(n_actions);
	action_lower_bound.Resize(n_actions);
	action_upper_bound += 1;
	
	state_action_lower_bound.Resize(n_states + n_actions);
	state_action_upper_bound.Resize(n_states + n_actions);
    
    reward = -1;
	
    endsim = false;
	Reset();
}

Acrobot::~Acrobot()
{
    // nothing to do
}

void Acrobot::Reset()
{
//	theta1		= urandom(state_lower_bound[0],state_upper_bound[0]);
//	theta2		= urandom(state_lower_bound[1],state_upper_bound[1]);
//	theta1Dot	= urandom(state_lower_bound[2],state_upper_bound[2]);
//	theta2Dot	= urandom(state_lower_bound[3],state_upper_bound[3]);
	theta1		= urandom() - 0.5;
	theta2		= urandom() - 0.5;
	theta1Dot	= urandom() - 0.5;
	theta2Dot	= urandom() - 0.5;

	state[0] = theta1;
	state[1] = theta2;
	state[2] = theta1Dot;
	state[3] = theta2Dot;
	endsim   = false;
	reward   = -1;
}

bool Acrobot::Act(const int& action)
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

void Acrobot::Simulate(const int action)
{
	real torque = action - 1.0;
	real d1;
	real d2;
	real phi_2;
	real phi_1;
	
	real theta2_ddot;
	real theta1_ddot;
	
	//torque is in [-1,1]
	//We'll make noise equal to at most +/- 1
	real theNoise = parameters.transitionNoise*2.0*(urandom() - 0.5);
	
	torque+=theNoise;
	
	int count = 0;
	while (!endsim && count < 4) {
		count++;
		
		d1 = parameters.m1 * pow(parameters.lc1, 2.0) + parameters.m2 * ( pow(parameters.l1, 2.0) + pow(parameters.lc2, 2.0) + 2.0 * parameters.l1 * parameters.lc2 * cos(theta2)) + parameters.I1 + parameters.I2;
		d2 = parameters.m2 * (pow(parameters.lc2, 2.0) + parameters.l1 * parameters.lc2 * cos(theta2)) + parameters.I2;
		
		phi_2 = parameters.m2 * parameters.lc2 * parameters.g * cos(theta1 + theta2 - M_PI / 2.0);
		phi_1 = -(parameters.m2 * parameters.l1 * parameters.lc2 * pow(theta2Dot, 2.0) * sin(theta2) - 2.0 * parameters.m2 * parameters.l1 * parameters.lc2 * theta1Dot * theta2Dot * sin(theta2)) + (parameters.m1 * parameters.lc1 + parameters.m2 * parameters.l1) * parameters.g * cos(theta1 - M_PI / 2.0) + phi_2;
		
		theta2_ddot = (torque + (d2 / d1) * phi_1 - parameters.m2 * parameters.l1 * parameters.lc2 * pow(theta1Dot, 2.0) * sin(theta2) - phi_2) / (parameters.m2 * pow(parameters.lc2, 2.0) + parameters.I2 - pow(d2, 2.0) / d1);
		theta1_ddot = -(d2 * theta2_ddot + phi_1) / d1;
		
		theta1Dot += theta1_ddot * parameters.dt;
		theta2Dot += theta2_ddot * parameters.dt;
		
		theta1 += theta1Dot * parameters.dt;
		theta2 += theta2Dot * parameters.dt;
	}
	if (abs(theta1Dot) > parameters.maxTheta1Dot) {
		theta1Dot = signum(theta1Dot) * parameters.maxTheta1Dot;
	}
	
	if (abs(theta2Dot) > parameters.maxTheta2Dot) {
		theta2Dot = signum(theta2Dot) * parameters.maxTheta2Dot;
	}
	/* Put a hard constraint on the Acrobot physics, thetas MUST be in [-PI,+PI]
	 * if they reach a top then angular velocity becomes zero
	 */
	if (abs(theta2) > M_PI) {
		theta2 = signum(theta2) * M_PI;
		theta2Dot = 0;
	}
	if (abs(theta1) > M_PI) {
		theta1 = signum(theta1) * M_PI;
		theta1Dot = 0;
	}

	real feet_height = -(parameters.l1 * cos(theta1) + parameters.l2 * cos(theta2));
	
	real firstJointEndHeight = parameters.l1 * cos(theta1);
	//Second Joint height (relative to first joint)
	real secondJointEndHeight = parameters.l2 * sin(M_PI / 2 - theta1 - theta2);
	feet_height = -(firstJointEndHeight + secondJointEndHeight);

	if (feet_height > parameters.acrobotGoalPosition) {
		endsim = true;
		reward = 0;
	}
	else {
		endsim = false;
		reward = -1;
	}
	state[0] = theta1;
	state[1] = theta2;
	state[2] = theta1Dot;
	state[3] = theta2Dot;
}

real Acrobot::signum(const real& num) {
	if(num == 0) {
		return 0.0;
	} else if(num > 0.0) {
		return 1.0;
	} else {
		return -1.0;
	}

}

