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

#include "UCRL2.h"
#include "ValueIteration.h"

/** Constructor.

    \param n_states_ number of states.
    \param n_actions_ number of actions.
    \param gamma_ discount factor.
    \param delta_ error probability to use in the bounds
    \rng 
*/
UCRL2::UCRL2(int n_states_,
             int n_actions_,
             real gamma_,
             DiscreteMDPCounts* model_,
             RandomNumberGenerator* rng_, 
			 real delta)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      confidence_interval(delta),
      model(model_),
	  rng(rng_),
      total_steps(0),
      update_interval(1),
      next_update(0),
      known_rewards(false),
      n_resets(0)
{
    state = -1;
    value_iteration = new OptimisticValueIteration(model, gamma);
    tmpQ.resize(n_actions);
}
UCRL2::~UCRL2()
{
    delete value_iteration;
}
void UCRL2::Reset()
{
    state = -1;
    n_resets++;
    //confidence_interval = 1.0 / (real) n_resets;
    //model->Reset();
}
/// Full observation
real UCRL2::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        model->AddTransition(state, action, reward, next_state);
    }
    return 0.0;
}
/// Partial observation 
real UCRL2::Observe (real reward, int next_state, int next_action)
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
int UCRL2::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);

    // update the model
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;
    total_steps++;
    if (total_steps >= next_update) {
        update_interval += n_states;
        next_update = total_steps + update_interval;
        //printf(" # next update: %d (interval %d)\n", next_update, update_interval);
        if (known_rewards) {
            value_iteration->ComputeStateValuesKnownRewards(confidence_interval, 1e-3, -1);
        } else {
            value_iteration->ComputeStateValues(confidence_interval, 1e-3, -1);
            //confidence_interval *= 0.5;
        }
        //const DiscreteMDP* mdp = model->getMeanMDP();
        //ValueIteration mean_vi(mdp, gamma);
        //mean_vi.ComputeStateValues(1e-6, -1);
        //printf("OVI\n");
        //value_iteration->getValues().print(stdout);
        //        printf("MVI\n");
        //mean_vi.getValues().print(stdout);
    }
    for (int i=0; i<n_actions; i++) {
        tmpQ[i] = value_iteration->getValue(next_state, i);
    }
    
    action = ArgMax(tmpQ);
    //printf("%f %d %d %f# rsaV'\n", reward, state, action, value_iteration->getValue(state, action));
    return action;
}


