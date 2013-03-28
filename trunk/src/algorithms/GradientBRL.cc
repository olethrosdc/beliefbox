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

#include "GradientBRL.h"

GradientBRL::GradientBRL(int n_states_,
                         int n_actions_,
                         real gamma_,
                         real epsilon_,
                         real step_size_,
                         MDPModel* model_,
                         RandomNumberGenerator* rng_,
                         bool use_upper_bound_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      epsilon(epsilon_),
      step_size(step_size_),
      model(model_),
      tmpQ(n_actions),
      rng(rng_),
      T(0),
      update_interval(n_states + 1),
      next_update(0),
      use_upper_bound(use_upper_bound_)
{
    printf("# Starting Gradient-B-RL with update interval %d\n",
           update_interval);

    current_state = -1;

    mdp_sample = model->generate();
    value_iteration = new ValueIteration(mdp_sample, gamma);
}

GradientBRL::~GradientBRL()
{
    delete value_iteration;
    delete mdp_sample;
}

void GradientBRL::Reset()
{
    current_state = -1;
    next_update = T;
    Resample();
}
/// Full observation
real GradientBRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        model->AddTransition(state, action, reward, next_state);
    }
    return 0.0;
}
/// Partial observation 
real GradientBRL::Observe (real reward, int next_state, int next_action)
{
    if (current_state >= 0) {
        model->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}


void GradientBRL::Resample()
{
    delete mdp_sample;
    mdp_sample = model->generate();
    value_iteration->setMDP(mdp_sample);
}



/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int GradientBRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);
    T++;

    // update the model
    if (current_state >= 0) {
        model->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;

    // Update MDPs
    //mdp_list[0] = model->getMeanMDP();
    // Do note waste much time generating MDPs
    
    bool do_update = false;
    if (T >= next_update) {    
        do_update = true;
    }
    if (do_update) {    
        update_interval += 1;
        next_update = T + update_interval;
        Resample();
    }

    real step_size_t = 1.0 / (1.0 / step_size + 0.001 * (real) T);

    if (use_upper_bound) {
        value_iteration->PartialUpdate(step_size_t);
    } else {
        value_iteration->PartialUpdateOnPolicy(step_size_t);
    }

    Vector tmpQ = value_iteration->getValues(current_state);

    int next_action;
    real epsilon_t = epsilon / (1.0 + sqrt((real) T));

    // choose action
    if (urandom()<epsilon_t) {
        next_action = rng->discrete_uniform(n_actions);
        //printf("RANDOM %d\n", next_action);
    } else {
        next_action = ArgMax(tmpQ);
    }
    current_action = next_action;
    //printf("%f %d #epsilon\n", epsilon_t, action);
    
    return current_action;
}


