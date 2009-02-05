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

#include "QLearning.h"

QLearning::QLearning(int n_states_,
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
      Q(n_states_, n_actions_),
      el(n_states_, n_actions_)
{
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = 0.0;//initial_value;
        }
    }
    exploration_policy->setValueMatrix(&Q);
    Reset();
}

void QLearning::Reset()
{
    state = -1;
    action = -1;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            el(s,a) = 0.0;
        }
    }

}

real QLearning::Observe (int action, int next_state, real reward)
{
    int a_max = 0;


    real Qa_max = Q(next_state, a_max);
    // select maximising action
    for (int i=1; i<n_actions; ++i) {
        if (Q(next_state, i) > Qa_max) {
            a_max = i;
            Qa_max = Q(next_state, a_max);
        }
    }

    real n_R = (reward - baseline) +  gamma*Qa_max; // partially observed return
    real TD = 0.0;
    if (state >= 0 && action >= 0) {
        real p_R = Q(state, action); // predicted return
        TD = n_R - p_R;
        
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_actions; ++j ) {
                el(i,j) *= lambda;
            }
        }
        el(state, action) = 1;
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_actions; ++j ) {
                Q(i, j) += el(i, j)*alpha*TD;	    
            }
        }
    }
    state = next_state; // fall back next state;
    return TD;
}

real QLearning::Observe (real reward, int next_state, int next_action)
{
    real TD = Observe(action, next_state, reward);
    action = next_action; // fall back next action
    return TD;
}

int QLearning::Act(real reward, int next_state)
{
    exploration_policy->Observe(reward, next_state);
    int next_action = exploration_policy->SelectAction();//Q, next_state);
    Observe(reward, next_state, next_action);
    return next_action;
}
