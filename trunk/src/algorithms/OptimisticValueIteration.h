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

#ifndef OPTIMISTIC_VALUE_ITERATION_H
#define OPTIMISTIC_VALUE_ITERATION_H

#include "DiscreteMDPCounts.h"
#include "real.h"
#include <vector>

class OptimisticValueIteration
{
public:
    DiscreteMDPCounts* mdp;
    real gamma;
    int n_states;
    int n_actions;
    std::vector<real> V;
    std::vector<real> dV;
    std::vector<real> pV;
    std::vector<real*> Q;
    std::vector<real> Q_data;
    std::vector<real*> dQ;
    std::vector<real> dQ_data;
    std::vector<real*> pQ;
    std::vector<real> pQ_data;
    real Delta;
    real baseline;
    OptimisticValueIteration(DiscreteMDPCounts* mdp, real gamma, real baseline=0.0);
    ~OptimisticValueIteration();
    void Reset();
    void ComputeStateValues(real epsilon, real threshold, int max_iter=-1);
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

