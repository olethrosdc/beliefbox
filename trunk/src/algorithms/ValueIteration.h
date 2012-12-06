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

#ifndef VALUE_ITERATION_H
#define VALUE_ITERATION_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "Vector.h"
#include "real.h"
#include <vector>

/** A value iteration algorithm for discrete MDPs */
class ValueIteration
{
protected:
    const DiscreteMDP* mdp; ///< pointer to the MDP
public:
    real gamma; ///< discount factor
    int n_states; ///< number of states
    int n_actions; ///< number of actions
    Vector V; ///< state values
    Vector dV; ///< difference between state values
    Vector pV; ///< previous statate value
    Matrix Q; ///< state-action value
    Matrix dQ; ///< state-action value difference
    Matrix pQ; ///< previous state-action values
    real Delta;
    real baseline;
    ValueIteration(const DiscreteMDP* mdp, real gamma, real baseline=0.0);
    ~ValueIteration();
    void Reset();
    inline void ComputeStateValues(real threshold, int max_iter=-1)
    {
        ComputeStateValuesElimination(threshold, max_iter);
        //ComputeStateValuesStandard(threshold, max_iter);
    }
    void ComputeStateValuesStandard(real threshold, int max_iter=-1);
    void ComputeStateValuesAsynchronous(real threshold, int max_iter=-1);
    void ComputeStateValuesElimination(real threshold, int max_iter=-1);
    void ComputeStateActionValues(real threshold, int max_iter=-1);
    /// Set the MDP to something else
    inline void setMDP(const DiscreteMDP* mdp_)
    {
        mdp = mdp_;
    }
    inline real getValue (int state, int action)
    {
        return Q(state, action);
    }
    inline real getValue (int state)
    {
        return V(state);
    }
    FixedDiscretePolicy* getPolicy();
    inline Matrix getValues() const
    {
        return Q;
    }
    inline Vector getStateValues() const
    {
        return V;
    }
    inline Vector getValues(int s) const
    {
        assert(s >= 0 && s < n_states);
        return Q.getRow(s);
    }
};
#endif

