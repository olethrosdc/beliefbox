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
                     real initial_value_,
                     real baseline_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      lambda(lambda_),
      alpha(alpha_),
      initial_value(initial_value_),
      baseline(baseline_),
      Q(n_states_, n_actions_),
      el(n_states_, n_actions_)
{
    state = -1;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = initial_value;
            el(s,a) = 0.0;
        }
    }
}

real QLearning::observe (int a, int ns, real r)
{
    real p_r = Q(state, action); // predicted reward
    real TD = (r - baseline) - p_r; // prediction error

    Q(s, a) = p_r + alpha * TD;
    
}
