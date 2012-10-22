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

#ifndef PUDDLEWORLD_H
#define PUDDLEWORLD_H

#include "Environment.h"
#include "Matrix.h"
#include "Vector.h"
#include "real.h"

/**The Puddle world environment.*/
class PuddleWorld: public Environment<Vector, int>
{
protected:
	struct Parameters{
		real	U_POS_X;		///< Upper bound on axis x
		real	L_POS_X;		///< Lower bound on axis x
		real	U_POS_Y;		///< Upper bound on axis y
		real	L_POS_Y;		///< Lower bound on axis y
		int 	NUMPUDDLES;		///< Number of puddles
		Matrix	U_POS_P;		///< Upper point on puddles
		Matrix	L_POS_P;		///< Lower point on puddles
		Vector	RADIUSPUDDLES;	///< Radius of puddles
		real	AGENTSPEED;		///< Step length
		real	MCNOISE;		///< Transition noise
	};
	Parameters default_parameters;
	Parameters parameters;
	Vector state_action_upper_bound;
	Vector state_action_lower_bound;
	Vector action_upper_bound;
	Vector action_lower_bound;
	Parameters SetDefaultParameters(){
		Parameters p;
		Matrix U_POS_P(2,2);
		Matrix L_POS_P(2,2);
		Vector RADIUSPUDDLES(2);
		p.U_POS_X		= 1.0;
		p.L_POS_X		= 0.0;
		p.U_POS_Y		= 1.0;
		p.L_POS_Y		= 0.0;
		p.NUMPUDDLES	= 2;
		p.U_POS_P		= U_POS_P;
		p.U_POS_P(0,0)	= 0.1;
		p.U_POS_P(0,1)	= 0.75;
		p.U_POS_P(1,0)	= 0.45;
		p.U_POS_P(1,1)	= 0.4;
		p.L_POS_P		= L_POS_P;
		p.L_POS_P(0,0)	= 0.45;
		p.L_POS_P(0,1)	= 0.75;
		p.L_POS_P(1,0)	= 0.45;
		p.L_POS_P(1,1)	= 0.8;
		p.RADIUSPUDDLES = RADIUSPUDDLES;
		p.RADIUSPUDDLES[0]	= 0.1;
		p.RADIUSPUDDLES[1]	= 0.1;
		p.AGENTSPEED		= 0.05;
		p.MCNOISE			= 0.1;
		return p;
	};

	void Simulate();
	bool endsim;
public:
	PuddleWorld(bool random_parameters = false);
	virtual ~PuddleWorld();
	virtual void Reset();
	virtual bool Act(const int action);
	virtual void Simulate(const int action);
	real DistPointToPuddle(const int puddle);

	
	Vector& StateActionUpperBound(){
		return state_action_upper_bound;
	}
	Vector& StateActionLowerBound(){
		return state_action_lower_bound;
	}
	Vector& ActionUpperBound(){
		return action_upper_bound;
	}
	Vector& ActionLowerBound(){
		return action_lower_bound;
	}
	
	virtual void setRandomness(real randomness){
		parameters.MCNOISE = randomness;
	}
};

#endif

	
