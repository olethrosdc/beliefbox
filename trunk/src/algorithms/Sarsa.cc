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

#include "Sarsa.h"

Sarsa::Sarsa(int n_states_,
             int n_actions_,
             real gamma_,
             real lambda_,
             real alpha_,
             VFExplorationPolicy* exploration_policy_,
             real initial_value_,
             real baseline_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      lambda(lambda_),
      alpha(alpha_),
      exploration_policy(exploration_policy_),
      initial_value(initial_value_),
      baseline(baseline_),
      Q(n_states_, n_actions_, Matrix::CHECK_BOUNDS),
      el(n_states_, n_actions_, Matrix::CHECK_BOUNDS)
{
    assert (lambda >= 0 && lambda <= 1);
    assert (alpha >= 0 && alpha <= 1);
    assert (gamma >=0 && gamma <= 1);
    
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = initial_value;
        }
    }
    exploration_policy->setValueMatrix(&Q);
    Reset();
}


Sarsa::~Sarsa() 
{
}


void Sarsa::Reset()
{
    state = -1;
    action = -1;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            el(s,a) = 0.0;
        }
    }

}

real Sarsa::Observe (int state, int action, real reward, int next_state, int next_action)
{
    real n_R = (reward - baseline) + gamma*Q(next_state, next_action); // partially observed return
    real p_R = Q(state, action); // predicted return
    real TD = n_R - p_R;

    Q(state, action) += alpha * TD;
    
    return TD;
}



real Sarsa::Observe (real reward, int next_state, int next_action)
{
    real n_R = (reward - baseline) + gamma*Q(next_state, next_action); // partially observed return
    //real p_R = 0.0;
    real TD = 0.0;
    if (state >= 0 && action >= 0) {
        real p_R = Q(state, action); // predicted return
        TD = n_R - p_R;
    

        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_actions; ++j ) {
                el(i,j) *= lambda;
            }
        }

        if (state >= 0 && action >= 0) {
            el(state, action) = 1;
        }

        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_actions; ++j ) {
                Q(i, j) += el(i, j) * alpha * TD;
            }
        }
    }
    state = next_state; // fall back next state;
    action = next_action; // fall back next action
    
    return TD;
}

int Sarsa::Act(real reward, int next_state)
{
    exploration_policy->Observe(reward, next_state);
    exploration_policy->setValueMatrix(&Q);
    int next_action = exploration_policy->SelectAction();
    Observe(reward, next_state, next_action);
    //printf ("Sarsa: %f %d %d\n", reward, next_state, next_action);
    //Q.print(stdout);
    return next_action;
}
