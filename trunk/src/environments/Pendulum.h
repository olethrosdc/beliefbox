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

#ifndef PENDULUM_H
#define PENDULUM_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

class Pendulum : public Environment<Vector, int>
{
protected:
    const real pendulum_mass; 		///< pendulum mass (kg)
    const real cart_mass;              ///< cart mass (kg)
    const real pendulum_length;        ///< pendulum length (m)
    const real gravity;  	        ///< gravity constant (g)
    const real max_noise;                  ///< noise
    const real Dt;                     ///< time constant
    const real CCa;		       ///< inverse total mass  
    static const int n_states = 2;     // state dimensions
    static const int n_actions = 3;     // action dimensions
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
    void Simulate();
    bool endsim;
    void penddot(Vector& xdot, real u, Vector& x);
    void pendulum_simulate(int action);
public:
    Pendulum(real pendulum_mass_ = 2.0,
             real cart_mass_ = 8.0,
             real pendulum_length_ = 0.5,
             real gravity_ = 9.8,
             real max_noise_ = 0.01,//10.0,
             real Dt_ = 0.01);
    virtual ~Pendulum();
    virtual void Reset();
    virtual bool Act(const int action);
    virtual void Simulate(const int action);
    Vector& StateUpperBound()
    {
        return state_upper_bound;
    }
    Vector& StateLowerBound()
    {
        return state_lower_bound;
    }
    Vector& StateActionUpperBound()
    {
        return state_action_upper_bound;
    }
    Vector& StateActionLowerBound()
    {
        return state_action_lower_bound;
    }
    Vector& ActionUpperBound()
    {
        return action_upper_bound;
    }
    Vector& ActionLowerBound()
    {
        return action_lower_bound;
    }

};








#endif
