// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_POMDP_H
#define DISCRETE_POMDP_H

#include "Matrix.h"


class DiscretePOMDP
{
protected:
    int n_states;
    int n_obs;
    int n_actions;
    Matrix Transitions;
    Matrix Observations;
    int state;
    int observation;
    real reward;
public:
    DiscretePOMDP(int n_states_, int n_obs_, int n_actions_);

    // use these to set the current state of the POMDP
    void setObservation(int x)
    {
        observation = x;
    }
    void setReward(real r)
    {
        reward = r;
    }
    void setState(int s)
    {
        state = s;
    }

    // get information about the POMDP
    int getNStates()
    {
        return n_states;
    }
    int getNActions()
    {
        return n_actions;
    }
    int getNObservations()
    {
        return n_obs;
    }
    int getObservation()
    {
        return observation;
    }
    real getNextStateProbability(int state, int action, int next_state) const
    {
        return Transitions(state*n_actions + action, next_state);
    }
	real getObservationProbability(int state, int action, int observation) const
    {
        return Observations(state*n_actions + action, observation);
    }

    // change the transition/observation matrices
    void setNextStateProbability(int state, int action, int next_state, real p)
    {
        Transitions(state*n_actions + action, next_state) = p;
    }
	void setObservationProbability(int state, int action, int observation, real p)
    {
        Observations(state*n_actions + action, observation) = p;
    }

    /// Check that the POMDP parameters are sane
    void check();
};

#endif
