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

#include "SarsaDirichlet.h"
#include <limits>

SarsaDirichlet::SarsaDirichlet(int n_states_,
                               int n_actions_,
                               real gamma_,
                               real lambda_,
                               real alpha_,
                               VFExplorationPolicy* exploration_policy_,
                               real initial_value_,
                               real baseline_)
    : Sarsa(n_states_,
            n_actions_,
            gamma_,
            lambda_,
            alpha_,
            exploration_policy_,
            initial_value_,
            baseline_),
      C(n_states_ * n_actions_, n_states)
{
    real lambda_count;
    if (lambda == 0.0) {
        Swarning("Lambda is 0\n");
        lambda_count = (real) std::numeric_limits<int>::max();
    } else { 
        lambda_count = (real) n_states / (lambda) - (real) n_states;
    }

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int i = s*n_actions + a;
            for (int s2=0; s<n_states; s++) {
                C(i, s2) = lambda_count;
            }
        }
    }
    exploration_policy->setValueMatrix(&Q);
    Reset();
}


SarsaDirichlet::~SarsaDirichlet() 
{
}


real SarsaDirichlet::Observe (real reward, int next_state, int next_action)
{
    real n_R = (reward - baseline) + gamma*Q(next_state, next_action); // partially observed return
    real p_R = 0.0;
    real TD = 0.0;
    if (state >= 0 && action >= 0) {
        p_R = Q(state, action); // predicted return
        TD = n_R - p_R;

        // find lambda
        // model assumes P(r|s',s,a) = P(r|s,a)
        real lambda_t = 0.0; // default value for initial state thing
        if (state >= 0) {
            int i = state*n_actions + action;
            real Z = (C.RowSum(i));
            lambda_t = C(i, next_state);
            C(i, next_state) += 1.0;
            if (Z > 0) {
                lambda_t /= Z;
            }
        }

        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_actions; ++j ) {
                el(i,j) *= lambda_t;
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

int SarsaDirichlet::Act(real reward, int next_state)
{
    exploration_policy->Observe(reward, next_state);
    int next_action = exploration_policy->SelectAction();//Q, next_state);
    Observe(reward, next_state, next_action);
    return next_action;
}
