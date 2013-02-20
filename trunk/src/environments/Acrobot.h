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

#ifndef ACROBOT_H
#define ACROBOT_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"
#include "AbstractPolicy.h"
#include "Random.h"

class Acrobot:public Environment<Vector, int>
{
protected:
    struct Parameters {
        real maxTheta1; 		
        real maxTheta2;           
        real maxTheta1Dot;       
        real maxTheta2Dot;  	            ///< gravity constant (g)
        real m1;             
        real m2;                    ///< time constant
		real l1;   
		real l2;
		real lc1;
		real lc2;
		real I1;
		real I2;
		real g;
		real dt;
		real acrobotGoalPosition;
		real transitionNoise;   ///tyre radious
    };
	static Parameters default_parameters;
    Parameters parameters;
	
	static const int n_states = 4;     // state dimensions
    static const int n_actions = 3;     // action dimensions
    real theta1, theta2, theta1Dot, theta2Dot;
	Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
	real signum(const real& num);
public:
	Acrobot(bool random_parameters = false);
	virtual ~Acrobot();
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
        parameters.transitionNoise = randomness;
    }	
	virtual const char* Name() const
    {
        return "Acrobot";
    }	
	void Show()
	{
		printf("%f %f %f %f %f %f # params (Bike)\n",
			   parameters.m1,
			   parameters.m2,
			   parameters.l1,
			   parameters.l2,
			   parameters.lc1,
			   parameters.lc2);
	}
	
};

class AcrobotGenerator
	{
	public:
		Acrobot Generate(bool random = true)
		{
			return Acrobot(random);
		}
	};

#endif
