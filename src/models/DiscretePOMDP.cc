// -*- Mode: c++ -*-
// copyright (c) 2009-2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscretePOMDP.h"

DiscretePOMDP::DiscretePOMDP(int n_states_, int n_obs_, int n_actions_)
    : n_states(n_states_),
      n_obs(n_obs_),
      n_actions(n_actions_),
      Transitions(n_states * n_actions, n_states),
      Observations(n_states * n_actions, n_obs)
{
    
    real p_state = 1.0 / (real) n_states;
    real p_obs = 1.0 / (real) n_obs;
    for (int i=0; i<n_states; ++i) { // uninitialised???
        for (int a=0; a<n_actions; ++a) {
            int k = i * n_actions + a;
            for (int j=0; j<n_states; ++j) { // uninitialised?
                Transitions(k, j) = p_state;
                Observations(i, j) = p_obs;
            }
        }
    }
}

void DiscretePOMDP::check()
{
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            real sum = 0;
            for (int s2=0; s2<n_states; ++s2)  {
                real p = getNextStateProbability(s, a, s2);
                sum += p;
                if (p < 0) {
                    Serror("P(s'=%d | s=%d, a=%d) = %f !\n", s2, s, a, p);
                    exit(-1);
                }
            }
            if (fabs(sum - 1) > 0.00001) {
                Serror("sum P(s' | s=%d, a=%d) = %f !\n", s, a, sum);
                exit(-1);
            }
        }
    }

}
