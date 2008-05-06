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
    assert (lambda >= 0 && lambda <= 1);
    assert (alpha >= 0 && alpha <= 1);
    assert (gamma >=0 && gamma <= 1);
    state = -1;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = initial_value;
            el(s,a) = 0.0;
        }
    }
}

int Sarsa::Observe (real r, int ns)
{
    int na = 0;

    // select next action
    real n_R = (r - baseline) + gamma*Q(ns, na); // partially observed return
    real p_R = Q(s,a); // predicted return
    real TD = n_R - p_R;

    for (int i=0; i<n_states; ++i) {
	for (int j=0; j<n_actions; ++j ) {
	    el(i,j) *= lambda;
	}
    }
    el(s, a) = 0;
    Q(s, a) = p_r + alpha * TD;
	 

    s = ns; // fall back next state;
}
