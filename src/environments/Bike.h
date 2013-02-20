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

#ifndef BIKE_H
#define BIKE_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"
#include "AbstractPolicy.h"
#include "Random.h"

class Bike :public Environment<Vector, int>
{
protected:
    struct Parameters {
        real R1; 		
        real R2;           
        real R3;       
        real R_FACTOR;  	            ///< gravity constant (g)
        int	 N0_STATES2;             
        real dt;                    ///< time constant
		real v;   
		real g;
		real dCM;
		real c;
		real h;
		real Mc;
		real Md;
		real Mp;
		real M;
		real R;   ///tyre radious
		real l;   ///distance between the point where the front and back tyre touch the ground
		real pi;
		/// position of goal
		real x_goal;
		real y_goal;
		real radius_goal;
		real max_noise;   ///< max noise
    };
    static Parameters default_parameters;
    Parameters parameters;
	real sigma_dot;
	real I_bike;
	real I_dc;
	real I_dv;
	real I_dl;
	real omega, omega_dot, omega_d_dot;
	real theta, theta_dot, theta_d_dot;
	real psi, psi_goal;
	real xf, yf, xb, yb;	///Tyre position
    static const int n_states = 6;     // state dimensions
    static const int n_actions = 5;     // action dimensions
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
    void Simulate();
	real sign(const real& num);
public:
	Bike(bool random_parameters = false);
	virtual ~Bike();
    virtual void Reset();
    virtual bool Act(const int& action);
    virtual void Simulate(const int action);	
    const Vector& StateUpperBound() const
    {
        return state_upper_bound;
    }
    const Vector& StateLowerBound() const
    {
        return state_lower_bound;
    }
    const Vector& StateActionUpperBound() const
    {
        return state_action_upper_bound;
    }
    const Vector& StateActionLowerBound() const
    {
        return state_action_lower_bound;
    }
    const Vector& ActionUpperBound() const
    {
        return action_upper_bound;
    }
    const Vector& ActionLowerBound() const
    {
        return action_lower_bound;
    }
	virtual real getTransitionProbability(const Vector& state, const int& action, const Vector& next_state) const
    {
        return 1.0; 
    }
	
    virtual real getExpectedReward(const Vector& state, const int& action) const 
    {
        return 0.0;
	}
	virtual void setRandomness(real randomness)
    {
        parameters.max_noise = randomness;
    }	
	virtual const char* Name() const
    {
        return "Bike";
    }	
	void Show()
	{
		printf("%f %f %f %f %d %f # params (Bike)\n",
			   parameters.R1,
			   parameters.R2,
			   parameters.R3,
			   parameters.R_FACTOR,
			   parameters.N0_STATES2,
			   parameters.dt);
	}
	float calc_dist_to_goal(float xf, float xb, float yf, float yb);
	float calc_angle_to_goal(float xf, float xb, float yf, float yb);
	
};

class BikeGenerator
{
public:
	Bike Generate()
	{
		return Bike(true);
	}
};

#endif
