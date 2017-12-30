// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "WeightedQLearning.h"

/** Initialise Weighted Q-learning.
	 
    The main parameters are n_models, which is the number of models to
    use, and epsilon, which is the probability with which each model
    is activated within a round. Always, only one model is used to
    make the decisions.

*/
WeightedQLearning::WeightedQLearning(int n_models_,
                                     int n_states_,
                                     int n_actions_,
                                     real gamma_,
                                     real alpha_,
                                     real epsilon_,
                                     real initial_value_,
                                     real baseline_)
    : n_models(n_models_),
      n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      alpha(alpha_),
      epsilon(epsilon_),
      initial_value(initial_value_),
      baseline(baseline_),
      current_model(0),
      update_interval(1),
      t(0),
      T(1)
{
    Q.resize(n_models);
    for (int m=0; m<n_models; m++) {
        Q[m].Resize(n_states, n_actions);
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++) {
                Q[m](s, a) = initial_value + urandom();
            }
        }
    }
    Reset();
}

/** Reset.
	
	Set the current state/action to invalid values. Clear eligibility traces.
*/
void WeightedQLearning::Reset()
{
    state = -1;
    action = -1;
    current_model = 0;
    update_interval = n_states;
    T = update_interval;
    t = 0;
}



/** Observe the current action and resulting next state and reward.

    We only need the next reward, state, and action pair, since the previous 
    state and action are saved by the algorithm.

	@param reward \f$r_{t+1}\f$	
	@param next_state \f$s_{t+1}\f$
	@param next_action \f$a_{t+1}\f$
*/
real WeightedQLearning::UpdateModel(int m, real reward, int next_state, int next_action)
{
    // select maximising action for the next state
    int a_max = 0;
    real Qa_max = Q[m](next_state, a_max);
    for (int i=1; i<n_actions; ++i) {
        if (Q[m](next_state, i) > Qa_max) {
            a_max = i;
            Qa_max = Q[m](next_state, a_max);
        }
    }

    real n_R = (reward - baseline) +  gamma*Qa_max; // partially observed return
    real p_R = Q[m](state, action); // predicted return
    real TD = n_R - p_R;
    real delta = alpha * TD;
    Q[m](state, action) += delta;	    
    //printf ("%d %f # m, TD\n", m, TD);
    return TD;
}

real WeightedQLearning::Observe(real reward, int next_state, int next_action)
{
    real P = 1;
    if (state >= 0 && action >= 0) {
        for (int i=0; i<n_models; ++i) {
            if (urandom() < P) {
                int m = urandom(0, n_models);
                UpdateModel(m, reward, next_state, next_action);
                P *= epsilon;
            } else {
                break;
            }
        }

    } 
    state = next_state; // fall back next state;
    action = next_action;
    return 0;
}

int WeightedQLearning::Act(real reward, int next_state)
{
    if (t > T) {
        current_model = urandom(0, n_models);
        update_interval += 1;
        T += update_interval;
    }
    int m = current_model;
    int next_action = 0;
    //printf ("%d\n", m);
    real Qmax = Q[m](next_state, next_action);
    for (int i=0; i<n_actions; ++i) {
        real Qsa = Q[m](next_state, i);
        if (Qsa > Qmax) {
            Qmax = Qsa;
            next_action = i;
        }
    }
    Observe(reward, next_state, next_action);
    action = next_action;
    //printf ("WeightedQLearning: %f %d %d\n", reward, next_state, next_action);
    //Q.print(stdout);
    return next_action;
}
