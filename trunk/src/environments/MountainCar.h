// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// copyright (c) 2003-2008 Michail G. Lagoudakis
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

class MountainCar : Environment<Vector, int>
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
    void Simulate();
    bool endsim;
public:
    MountainCar();
    virtual ~MountainCar();
    virtual void Reset();
    virtual bool Act(int action);
};

real MountainCar::U_POS = 0.5;                // Upper bound on position
real MountainCar::L_POS = -1.2;             // Lower bound on position
real MountainCar::U_VEL = 0.07;               // Upper bound on velocity
real MountainCar::L_VEL = -0.07;            // Lower bound on velocity
real MountainCar::INPUT = 0.001;           // contribution of input
real MountainCar::GRAVITY = 0.0025;          // contribution of gravity
real MountainCar::MCNOISE = 0.2






#endif
