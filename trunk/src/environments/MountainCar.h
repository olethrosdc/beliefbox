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
    static real U_POS; 		// Upper bound on position
    static real L_POS; 		// Lower bound on position
    static real U_VEL;  	        // Upper bound on velocity
    static real L_VEL;  		// Lower bound on velocity
    static real INPUT;		// contribution of input
    static real GRAVITY; 	// contribution of gravity
    static real MCNOISE;        // input noise
    static const int n_states = 2;     // state dimensions
    static const int n_actions = 3;     // action dimensions
    Vector state_upper_bound;
    Vector state_lower_bound;
    void Simulate();
    bool endsim;
public:
    MountainCar();
    virtual ~MountainCar();
    virtual void Reset();
    virtual bool Act(int action);
    virtual void Simulate(int action);
    Vector& StateUpperBound()
    {
        return state_upper_bound;
    }
    Vector& StateLowerBound()
    {
        return state_lower_bound;
    }
    void setNoise(real noise)
    {
        MCNOISE = noise;
    }
};








#endif
