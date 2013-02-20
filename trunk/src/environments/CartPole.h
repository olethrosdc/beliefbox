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

#ifndef CART_POLE_H
#define CART_POLE_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"
#include "AbstractPolicy.h"
#include "Random.h"

class CartPole : public Environment<Vector, int>
{
protected:
    struct Parameters {
        real GRAVITY; 		
        real MASSCART;         
        real MASSPOLE;       
        real LENGTH;             
        real FORCE_MAG;                   
		real TAU;
		real noise;
    };
    static Parameters default_parameters;
    Parameters parameters;
    static const real FOURTHIRDS		= 4.0/3.0;
    Vector state_action_upper_bound;
    Vector state_action_lower_bound;
    Vector action_upper_bound;
    Vector action_lower_bound;
	real TOTAL_MASS, POLEMASS_LENGTH;
    void Simulate();
    void penddot(Vector& xdot, real u, Vector& x);
    void pendulum_simulate(int action);
public:
    CartPole(bool random_parameters = false);
    virtual ~CartPole();
    virtual void Reset();
    virtual bool Act(const int& action);
    virtual void Simulate(const int action);
	virtual const char* Name() const
    {
        return "Cart Pole RL";
    }
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
        if (fabs(state[0]) > 0.2094) {
            return -1.0;
        } else {
            return 1.0;
        }
	}
	void Show()
	{
		printf("%f %f %f %f %f %f # params (CartPole)\n",
			   parameters.GRAVITY,
			   parameters.MASSCART,
			   parameters.MASSPOLE,
			   parameters.LENGTH,
			   parameters.FORCE_MAG,
			   parameters.TAU);
	}
	
};



class CartPoleGenerator
{
	public:
		CartPole Generate(bool random=true)
		{
			return CartPole(random);
		}
};



class HeuristicCartPolePolicy  : public AbstractPolicy<Vector, int>
{
protected:
	int n_actions;
    Vector state;
public:
	HeuristicCartPolePolicy()
	: n_actions(2), state(2)
	{
	}
	virtual ~HeuristicCartPolePolicy()
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
