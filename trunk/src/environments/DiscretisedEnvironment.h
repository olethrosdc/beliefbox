// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETISED_ENVIRONMENT_H
#define DISCRETISED_ENVIRONMENT_H

#include "Environment.h"
#include "Grid.h"

template<class T>
class DiscretisedEnvironment : public DiscreteEnvironment {
public:
    T& environment;
    int K;
    EvenGrid grid;
    DiscretisedEnvironment(T& environment_, int K_)
        : environment(environment_),
          K(K_),
          grid(environment.StateLowerBound(), environment.StateUpperBound(), K)
    {
        n_actions = environment.getNActions();
        n_states = grid.getNIntervals();
    }

    virtual ~DiscretisedEnvironment()
    {
    }

    virtual void Reset() {
        environment.Reset();
        state = grid.getInterval(environment.getState());
    }

    virtual bool Act(int action) {
        bool flag = environment.Act(action);
        state = grid.getInterval(environment.getState());
        reward = environment.getReward();
        return flag;
    }

    virtual const char* Name()
    {
        return "Discretised environment";
    }


};

#endif
