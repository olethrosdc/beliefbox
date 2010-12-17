// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SampleBasedRL.h"

SampleBasedRL::SampleBasedRL(int n_states_,
                           int n_actions_,
                           real gamma_,
                           real epsilon_,
                           MDPModel* model_,
                           int max_samples_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      epsilon(epsilon_),
      model(model_),
      max_samples(max_samples_),
      T(0),
      update_interval(10)
{
    state = -1;

    Vector w(max_samples);
    real w_i = 1.0 / (real) max_samples;
    mdp_list.resize(max_samples);
    mdp_list[0] = model->getMeanMDP();
    for (int i=1; i<max_samples; ++i) {
        mdp_list[i] = model->generate();
        w[i] = w_i;
    }

    value_iteration = new MultiMDPValueIteration(w, mdp_list, gamma);
    tmpQ.resize(n_actions);
}
SampleBasedRL::~SampleBasedRL()
{
    for (int i=1; i<max_samples; ++i) {
        delete mdp_list[i];
    }
    delete value_iteration;
}
void SampleBasedRL::Reset()
{
    state = -1;
    //model->Reset();
}
/// Full observation
real SampleBasedRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        model->AddTransition(state, action, reward, next_state);
    }
    return 0.0;
}
/// Partial observation 
real SampleBasedRL::Observe (real reward, int next_state, int next_action)
{
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;
    action = next_action;
    return 0.0;
}
/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int SampleBasedRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);
    T++;

    // update the model
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;

    // Update MDPs
    mdp_list[0] = model->getMeanMDP();
    // Do note waste much time generating MDPs
    if (T >= update_interval) {    
        update_interval = 2*T;
        for (int i=1; i<max_samples; ++i) {
            delete mdp_list[i];
            mdp_list[i] = model->generate();
            //model->ShowModel();
        }
    }
    
    // update values
    value_iteration->setMDPList(mdp_list);
    value_iteration->ComputeStateActionValues(0,1);
    for (int i=0; i<n_actions; i++) {
        tmpQ[i] = value_iteration->getValue(next_state, i);
        //printf ("Q[%d] = %f ", i, tmpQ[i]);
    }
    
    int next_action;
    // choose action
    if (urandom()<epsilon) {
        next_action = (int) floor(urandom(0.0, n_actions));
        //printf ("\n");
    } else {
        next_action = ArgMax(tmpQ);
    }
    action = next_action;

    return action;
}


