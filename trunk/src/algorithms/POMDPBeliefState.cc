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
DiscretePOMDPBeliefState::DiscretePOMDPBeliefState(DiscretePOMDP* pomdp_) : pomdp(pomdp_)
{
    n_states = pomdp->getNStates();
    belief.Resize(n_states);
    log_belief.Resize(n_states);
    Reset();
}

void DiscretePOMDPBeliefState::Reset() 
{
    real p = 1.0 / (real) n_states;
    real log_p = log(p);
    for (int i=0; i<n_states; ++i) {
        belief[i] = p;
        log_belief[i] = log_p;
    }
}

DiscretePOMDPBeliefState::~DiscretePOMDPBeliefState()
{
    
}




real DiscretePOMDPBeliefState::Observe(int a, int x, real r)
{
#if 1
    Vector next_belief(n_states);
    real sum = 0;
    for (int s=0; s<n_states; s++) {
        real p_obs = pomdp->getObservationProbability(s, a, x);
        next_belief[s] = belief[s] * p_obs;
        sum += next_belief[s];
    }
    for (int s=0; s<n_states; s++) {
        next_belief[s] /= sum;
    }

    // P(s'|x',a,b) = \sum_s P(s'|s,a) P(s|x',a,b) 
    for (int s2=0; s2<n_states; ++s2) {
        belief[s2] = 0;
        for (int s=0; s<n_states; ++s) {
            real p_ss = pomdp->getNextStateProbability(s, a, s2);
            belief[s2] += next_belief[s]* p_ss;
        }
    }
    for (int s2=0; s2<n_states; ++s2) {
        log_belief[s2] = log(belief[s2]);
    }


    return sum;

#else

    Vector log_belief2(n_states);
    real log_sum = LOG_ZERO;
    for (int s=0; s<n_states; s++) {
        real p_obs = pomdp->getObservationProbability(s, a, x);
        log_belief2[s] = log_belief[s] + log(p_obs);
        log_sum = logAdd(log_sum, log_belief2[s]);
    }
    for (int s=0; s<n_states; s++) {
        log_belief2[s] -= log_sum;
    }

    // P(s'|x',a,b) = \sum_s P(s'|s,a) P(s|x',a,b) 
    for (int s2=0; s2<n_states; ++s2) {
        log_belief[s2] = 0;
        for (int s=0; s<n_states; ++s) {
            real p_ss = pomdp->getNextStateProbability(s, a, s2);
            log_belief[s2] = logAdd(log_belief[s2], log_belief2[s] + log(p_ss));
        }
    }
    for (int s2=0; s2<n_states; ++s2) {
        belief[s2] = exp(log_belief[s2]);
        printf ("B(%d)=%f\n", s2, belief[s2]);
    }




    return exp(log_sum);


#endif
}

/** Calculate observation probability: \f$P(x_{t+1},r_{t+1} | b_t,a_t)\f$.

    First
    \f[
    P(s_{t+1} | a_t, b_t) = \sum_s P(s_{t+1}|a_t, s_t = s) P(s_t = s | b_t)
    \f]
    \f[
    P(x_{t+1} | a_t, b_t) = \sum_s P(x_{t+1} | s_{t+1} = s) P(s_{t+1} = s | a_t, b_t)
    \f]
    
 */
real DiscretePOMDPBeliefState::ObservationProbability(int a, int x, real r)
{
    real sum = 0;
    for (int s=0; s<n_states; s++) {
        real p_obs = pomdp->getObservationProbability(s, a, x);
        sum +=  belief[s] * p_obs;
    }
    return sum;
}


