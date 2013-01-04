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
#include "AbstractPolicy.h"
#include "Random.h"

class Pendulum : public Environment<Vector, int>
{
protected:
    struct Parameters {
        real pendulum_mass; 		///< pendulum mass (kg)
        real cart_mass;             ///< cart mass (kg)
        real pendulum_length;       ///< pendulum length (m)
        real gravity;  	            ///< gravity constant (g)
        real max_noise;             ///< noise
        real Dt;                    ///< time constant
    };
    static Parameters default_parameters;
    Parameters parameters;
    real CCa;		            ///< inverse total mass  
    static const int n_states = 2;     // state dimensions
    static const int n_actions = 3;     // action dimensions
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
    void Simulate();
    void penddot(Vector& xdot, real u, Vector& x);
    void pendulum_simulate(int action);
public:
    Pendulum(bool random_parameters = false);
    virtual ~Pendulum();
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
        if (fabs(state[0]) > M_PI/2.0) {
            return -1.0;
        } else {
            return 0.0;
        }
	}
	void Show()
	{
		printf("%f %f %f %f %f %f # params (Pendulum)\n",
			   parameters.pendulum_mass,
			   parameters.cart_mass,
			   parameters.pendulum_length,
			   parameters.gravity,
			   parameters.max_noise,
			   parameters.Dt);
	}

};



class PendulumGenerator
{
public:
    Pendulum Generate()
    {
        return Pendulum(true);
    }
};



class HeuristicPendulumPolicy  : public AbstractPolicy<Vector, int>
{
protected:
	int n_actions;
    Vector state;
public:
	HeuristicPendulumPolicy()
		: n_actions(2),
          state(2)
	{
	}
	virtual ~HeuristicPendulumPolicy()
	{
	}
	virtual int SelectAction()
	{
        int action = 1;
        if (urandom() < 0.1) {
            return rand()%3;
        }
        if (state[0] + state[1] > 0) {
            action = 2;
        }
        if (state[0] + state[1] < 0) {
            action = 0;
        }
        return action;
	}
	virtual void Observe (const Vector& previous_state, const int& action, real r, const Vector& next_state) 
	{
        state = next_state;
	}
    virtual void Observe (real r, const Vector& next_state) 
	{
        state = next_state;
	}
	virtual void Reset() 
	{
        state(0) = 0;
        state(1) = 0;
	}
};

#endif
