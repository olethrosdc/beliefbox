// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MOUNTAINCAR_3D_H
#define MOUNTAINCAR_3D_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

/** The mountain car environment.
 */
class MountainCar3D : public Environment<Vector, int>
{
protected:
    static real U_POS; 		///< Upper bound on position
    static real L_POS; 		///< Lower bound on position
    static real U_VEL;  	        ///< Upper bound on velocity
    static real L_VEL;  		///< Lower bound on velocity
    static real INPUT;		///< contribution of input
    static real GRAVITY; 	///< contribution of gravity
    static real MCNOISE;        ///< input noise
    Vector state_upper_bound;
    Vector state_lower_bound;
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
    void Simulate();
public:
    MountainCar3D();
    virtual ~MountainCar3D();
    virtual void Reset();
    virtual bool Act(int action);
    virtual void Simulate(int action);
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

    void setNoise(real noise)
    {
        MCNOISE = noise;
    }
};








#endif
