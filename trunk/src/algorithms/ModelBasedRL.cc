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

#include "ModelBasedRL.h"

ModelBasedRL::ModelBasedRL(int n_states_,
                           int n_actions_,
                           real gamma_,
                           real epsilon_,
                           MDPModel* model_,
						   RandomNumberGenerator* rng_,
                           bool use_value_iteration_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      epsilon(epsilon_),
      model(model_),
	  rng(rng_),
      use_value_iteration(use_value_iteration_),
      total_steps(0)
{
    state = -1;
    mdp = model->CreateMDP();
	assert(mdp);
    value_iteration = new ValueIteration(mdp, gamma);
    tmpQ.resize(n_actions);
}
ModelBasedRL::~ModelBasedRL()
{
    //mdp->ShowModel();
	delete mdp;
    delete value_iteration;
}
void ModelBasedRL::Reset()
{
    state = -1;
    //model->Reset();
}
/// Full observation
real ModelBasedRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        model->AddTransition(state, action, reward, next_state);
    }
    return 0.0;
}
/// Partial observation 
real ModelBasedRL::Observe (real reward, int next_state, int next_action)
{
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;
    action = next_action;
    total_steps++;
    return 0.0;
}
/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int ModelBasedRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);

    // update the model
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;

    if (use_value_iteration) {
		if (mdp) {
			delete mdp;
		}
        mdp = model->getMeanMDP();
        
        // update values
        value_iteration->setMDP(mdp);
        value_iteration->ComputeStateValues(1e-6,1);
        for (int i=0; i<n_actions; i++) {
            tmpQ[i] = value_iteration->getValue(next_state, i);
        }
    } else {
        exit(-1);
        for (int i=0; i<n_actions; i++) {
            tmpQ[i] = model->getExpectedReward(state, i);
        } 
    }
    
    // choose action
    int next_action;
    real epsilon_t = epsilon / (1.0+ sqrt((real) total_steps));
    total_steps++;
    if (rng->uniform()<epsilon_t) {
        next_action = rng->discrete_uniform(n_actions);
        //printf ("\n");
        //printf ("r ");
    } else {
        next_action = ArgMax(tmpQ);
        //printf ("m ");
        //printf (" * A = %d\n", next_action);
    }
    action = next_action;

    //printf("%d %f %d #epsilon\n", state, epsilon_t, action);
    return action;
}


