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
             MDPModel* model_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      model(model_)
{
    mdp = model->CreateMDP();
    value_iteration = new ValueIteration(mdp, gamma);
    tmpQ.resize(n_actions);
}
ModelBasedRL::~ModelBasedRL()
{
    delete value_iteration;
    delete mdp;
}
void ModelBasedRL::Reset()
{
    model->Reset();
}
/// Full observation
real ModelBasedRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    model->AddTransition(state, action, reward, next_state);
    return 0.0;
}
/// Partial observation 
real ModelBasedRL::Observe (real reward, int next_state, int next_action)
{
    model->AddTransition(state, action, reward, next_state);
    state = next_state;
    return 0.0;
}
/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int ModelBasedRL::Act(real reward, int next_state)
{

    DiscreteMDP* mdp = model->CreateMDP();
    value_iteration->mdp = mdp;
    value_iteration->ComputeStateActionValues(0.00, 1);
    for (int i=0; i<n_actions; i++) {
        tmpQ[i] = getValue(state, i);
    }
    int next_action = ArgMax(tmpQ);
    Observe(reward, next_state, action);
    state = next_state;
    action = next_action;
    return action;
}


