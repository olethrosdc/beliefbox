// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef AVERAGE_VALUE_ITERATION_H
#define AVERAGE_VALUE_ITERATION_H

#include "ValueIteration.h"
#include "DiscreteMDP.h"
#include "real.h"
#include <vector>

class AverageValueIteration
{
public:
    const DiscreteMDP* mdp;
    real gamma;
    int n_states;
    int n_actions;
    std::vector<real> V; ///< current value
    std::vector<real> dV; ///< value difference
    std::vector<real> pV; ///< previous value
    std::vector<real*> Q; ///< current action-state value
    std::vector<real> Q_data; ///< data
    std::vector<real*> dQ; ///< current action-state value difference
    std::vector<real> dQ_data; ///< current action-state value data
    std::vector<real*> pQ; ///< previous action-state value
    std::vector<real> pQ_data; ///< previous action-state value data
    std::vector<real> p_b; ///< stationary distribution for baseline
    real Delta; ///< difference norm for termination
    real baseline; ///< baseline for fixing overflow
    AverageValueIteration(const DiscreteMDP* mdp, real baseline=0.0);
    ~AverageValueIteration();
    void Reset();
    void ComputeStateValues(real threshold, int max_iter=-1);
    void ComputeStateActionValues(real threshold, int max_iter=-1);
    inline real getValue (int state, int action)
    {
        assert(state>=0 && state < n_states);
        assert(action>=0 && action < n_actions);
        return Q[state][action];
    }
    inline real getValue (int state)
    {
        assert(state>=0 && state < n_states);
        return V[state];
    }
    
};

#endif

