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

#ifndef POMDP_BELIEF_STATE_H
#define POMDP_BELIEF_STATE_H

#include "DiscretePOMDP.h"

class DiscreretePOMDPBeliefState
{
protected:
    DiscretePOMDP* pomdp;
    int n_states;
    Vector belief;
    Vector log_belief;
public: 
    DiscreretePOMDPBeliefState(DiscretePOMDP* pomdp_);
    ~DiscreretePOMDPBeliefState();

    real Observe(int a, int x, real r);

};

#endif
