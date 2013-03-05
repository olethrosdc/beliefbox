/*
 *  LinearDynamicQuadratic.h
 *  
 *
 *  Created by Poseidon on 20/12/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LINEARDYNAMICQUADRATIC_H
#define LINEARDYNAMICQUADRATIC_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

class LinearDynamicQuadratic:public Environment<Vector, int>
{
protected:
	struct Parameters {
        real U_POS;         ///< Upper bound on position
        real L_POS;         ///< Lower bound on position
        real U_VEL;         ///< Upper bound on velocity
        real L_VEL;         ///< Lower bound on velocity
        real INPUT;         ///< contribution of input
        real MCNOISE;       ///< input noise        
    };
    static Parameters default_parameters;
    Parameters parameters;
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
public:
	LinearDynamicQuadratic(bool random_parameters = false);
	virtual ~LinearDynamicQuadratic();
	virtual void Reset();
	virtual bool Act(const int& action);
	virtual void Simulate(const int& action);
	
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
	
	virtual void setRandomness( real randomness )
	{
		parameters.MCNOISE = randomness;
	}
	
	virtual const char* Name() const
	{
		return "LinearDynamicQuadratic";
	}
	
	void Show()
	{
		//printf("%f %f %f %f %f %f %f # params (LinearDynamicQuadratic)\n",
//			   parameters.U_POS,
//			   parameters.L_POS,
//			   parameters.U_VEL,
//			   parameters.L_VEL,
//			   parameters.INPUT,
//			   parameters.MCNOISE);
	}
};

#endif

