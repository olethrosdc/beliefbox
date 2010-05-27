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

#include "Pendulum.h"
#include "Random.h"

Pendulum::Pendulum(real pendulum_mass_,
                   real cart_mass_,
                   real pendulum_length_,
                   real gravity_,
                   real max_noise_,
                   real Dt_)
    : Environment<Vector, int>(2, 3),
      pendulum_mass(pendulum_mass_),
      cart_mass(cart_mass_),
      pendulum_length(pendulum_length_),
      gravity(gravity_), 
      max_noise(max_noise_),
      Dt(Dt_),
      CCa (1.0 / (pendulum_mass + cart_mass))
{
    assert(pendulum_mass > 0);
    assert(cart_mass > 0);
    assert(pendulum_length > 0);
    assert(gravity > 0);
    assert(max_noise > 0);
    assert(Dt > 0);
    state.Resize(2);
    state_upper_bound.Resize(2);
    state_lower_bound.Resize(2);
    state_upper_bound[0] = 4;
    state_upper_bound[1] = 10;
    state_lower_bound[0] = -4;
    state_lower_bound[1] = -10;
    endsim = false;
}

Pendulum::~Pendulum()
{
    // nothing to do
}

void Pendulum::Reset()
{

    // Theta
    state[0] =  0.0;//urandom(-0.01, 0.01);

    // dTheta/dt
    state[1] = 0.0;//urandom(-0.001, 0.001);
    endsim = false;
}

void Pendulum::penddot(Vector& xdot, real u, Vector& x)
{
    // Nonlinear model 
     
    double cx = cos(x[0]);
    real dtheta2 =x[1]*x[1];
    xdot[0] = x[1]; 
    xdot[1] = (gravity * sin(x[0]) - 
               0.5*CCa * pendulum_mass * pendulum_length * dtheta2 * sin(2.0*x[0]) -  CCa * cos(x[0]) * u ) / 
        ( 4.0/3.0*pendulum_length - CCa*pendulum_mass*pendulum_length*cx*cx ); 
   
}

bool Pendulum::Act(int action)
{
    // make sure we tell the guy we have terminated
    if (endsim) {
        reward = 0.0;
        return false;
    }
    
    // run
    Simulate(action);
    return true;
}

void Pendulum::Simulate(int action)
{
    Vector xdot(2);
    real input=0.0, noise, t;

  
    switch(action) {
    case 0: input = -50.0; break;
    case 1: input = 0.0; break;
    case 2:  input = +50.0; break;
    }

    noise = urandom(-max_noise, max_noise);
    input += noise;

    // Simulate for 0.1 seconds
    for (t=0.0; t<=0.1; t+=Dt) {

        penddot(xdot, input, state);
    
        state[0] += xdot[0] * Dt;
        state[1] += xdot[1] * Dt;

    }
  
    if (fabs(state[0]) > M_PI/2.0) {
        reward = -1.0;
        endsim = true;
    } else {
        reward = 0.0;
        endsim = false;
    }
  
   
}
