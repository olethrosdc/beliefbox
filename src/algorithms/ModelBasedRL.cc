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
                           MDPModel* model_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      epsilon(epsilon_),
      model(model_)
{
    state = -1;
    mdp = model->CreateMDP();
    value_iteration = new ValueIteration(mdp, gamma);
    tmpQ.resize(n_actions);
}
ModelBasedRL::~ModelBasedRL()
{
    delete value_iteration;
    if (mdp) {
        delete mdp;
    }
}
void ModelBasedRL::Reset()
{
    state = -1;
    model->Reset();
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
    return 0.0;
}
/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int ModelBasedRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);

    if (mdp) {
        delete mdp;
        mdp = NULL;
    }
    mdp = model->CreateMDP();
    //mdp->Check();
#if 0
    if (value_iteration) {
        delete value_iteration;
    }
    value_iteration = new ValueIteration(mdp, gamma);
    //value_iteration->ComputeStateValues(0.001, 1000000);
    value_iteration->ComputeStateActionValues(0.001, 1000000);
#else
    value_iteration->mdp = mdp;
    //value_iteration->ComputeStateValues(0.00, 1);
    value_iteration->ComputeStateActionValues(0.00, 1);
#endif
    //printf ("V(%d) = %f\n", next_state, value_iteration->getValue(next_state));
    for (int i=0; i<n_actions; i++) {
        tmpQ[i] = value_iteration->getValue(next_state, i);
        //printf ("Q(%d %d) = %f\n", next_state, i, tmpQ[i]);
    }
    int next_action;
    if (urandom()<epsilon) {
        next_action = (int) floor(urandom(0.0, n_actions));
    } else {
        next_action = ArgMax(tmpQ);
    }
    Observe(reward, next_state, action);
    state = next_state;
    action = next_action;
    return action;
}


