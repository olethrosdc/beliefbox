// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GOAL_STATE_ITERATION_H
#define GOAL_STATE_ITERATION_H

#include "DiscreteMDP.h"
#include "real.h"
#include <vector>


/** This algorithm calculates the distance between two states */
class GoalStateIteration
{
public:
    const DiscreteMDP* mdp;
    real gamma;
    int n_states;
    int n_actions;
    std::vector<real> V;
    std::vector<real> dV;
    std::vector<real> pV;
    real Delta;
    real baseline;
    GoalStateIteration(const DiscreteMDP* mdp);
    ~GoalStateIteration();
    void Reset();
    void ComputeStateValues(int goal_state, real threshold, int max_iter=-1);
    real GetMaximumDistanceFromState(int goal_state, real threshold, int max_iter=-1);
    real GetMDPDiameter(real threshold, int max_iter=-1);
    inline real getValue (int state)
    {
        assert(state>=0 && state < n_states);
        return V[state];
    }
    
};

#endif

