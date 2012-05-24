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

#include "DiscretePolicy.h"
#include "DiscreteMDPCounts.h"
#include "real.h"
#include <vector>



class OptimisticValueIteration
{
public:
    const DiscreteMDPCounts* mdp;
    real gamma;
    int n_states;
    int n_actions;
    Vector V; ///< state values
    Vector dV; ///< difference between state values
    Vector pV; ///< previous statate value
    Matrix Q; ///< state-action value
    Matrix dQ; ///< state-action value difference
    Matrix pQ; ///< previous state-action values
    real Delta;
    real baseline;

    OptimisticValueIteration(const DiscreteMDPCounts* mdp,
                             real gamma,
                             real baseline=0.0);
    ~OptimisticValueIteration();
    void Reset();
    void ComputeStateValues(real error_probability, real epsilon, real threshold, int max_iter=-1)
    {
        ComputeStateValuesStandard(error_probability,
                                   epsilon,
                                   threshold,
                                   max_iter);
    }
    void ComputeStateValuesStandard(real error_probability, real epsilon, real threshold, int max_iter=-1);
    inline real getValue (int state, int action)
    {
        return Q(state, action);
    }
    inline real getValue (int state)
    {
        return V(state);
    }
    FixedDiscretePolicy* getPolicy() const;
};

#endif

