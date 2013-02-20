// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// Adapted from code by Michail G. Lagoudakis, copyright (c) 2003-2008 
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MOUNTAINCAR_H
#define MOUNTAINCAR_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

/** The mountain car environment.
 */
class MountainCar : public Environment<Vector, int>
{
protected:
    struct Parameters {
        real U_POS;         ///< Upper bound on position
        real L_POS;         ///< Lower bound on position
        real U_VEL;         ///< Upper bound on velocity
        real L_VEL;         ///< Lower bound on velocity
        real INPUT;         ///< contribution of input
        real GRAVITY;       ///< contribution of gravity
        real MCNOISE;       ///< input noise        
    };
    static Parameters default_parameters;
    Parameters parameters;
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
    void Simulate();
	bool endsim;
public:
    MountainCar(bool random_parameters = false);
    virtual ~MountainCar();
    virtual void Reset();
    virtual bool Act(const int& action);
    virtual void Simulate(const int action);

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

    virtual void setRandomness(real randomness)
    {
        parameters.MCNOISE = randomness;
    }

	virtual const char* Name() const
    {
        return "Mountain Car";
    }

	void Show()
	{
		printf("%f %f %f %f %f %f %f # params (MountainCar)\n",
			   parameters.U_POS,
			   parameters.L_POS,
			   parameters.U_VEL,
			   parameters.L_VEL,
			   parameters.INPUT,
			   parameters.GRAVITY,
			   parameters.MCNOISE);
	}
    virtual real getTransitionProbability(const Vector& state, const int& action, const Vector& next_state) const
    {
        return 1.0;
    }

    virtual real getExpectedReward(const Vector& state, const int& action) const
    {
        return 0.0;
    }
};

class MountainCarGenerator
{
public:
    MountainCar Generate(bool random=true)
    {
      return MountainCar(random);
    }
};







#endif
