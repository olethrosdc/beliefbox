// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "POMDPBeliefState.h"

/** Construct a new belief state.
    
    Initialise the belief to uniform.
*/
DiscreretePOMDPBeliefState::DiscreretePOMDPBeliefState(DiscretePOMDP* pomdp_) : pomdp(pomdp_)
{
    n_states = pomdp->getNStates();
    belief.Resize(n_states);
    log_belief.Resize(n_states);
    real p = 1.0 / (real) n_states;
    real log_p = log(p);
    for (int i=0; i<n_states; ++i) {
        belief[i] = p;
        log_belief[i] = log_p;
    }
}

DiscreretePOMDPBeliefState::~DiscreretePOMDPBeliefState()
{
    
}


real DiscreretePOMDPBeliefState::Observe(int a, int x, real r)
{
    Vector next_belief(n_states);
    for (int s2=0; s2<n_states; ++s2) {
        next_belief[s2] = 0;
        for (int s=0; s<n_states; ++s) {
            next_belief[s2] += belief[s] * pomdp->getNextStateProbability(s, a, s2);
        }
    }
    real log_sum = LOG_ZERO;
    for (int s=0; s<n_states; s++) {
        log_belief[s] = log(next_belief[s]) + log(pomdp->getObservationProbability(s, x));
        log_sum = logAdd(log_sum, log_belief[s]);
    }

    for (int s=0; s<n_states; s++) {
        log_belief[s] -= log_sum;
        belief[s] = exp(log_belief[s]);
    }
    
    return exp(log_sum);
}

