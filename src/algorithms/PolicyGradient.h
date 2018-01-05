// -*- Mode: c++ -*-
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef POLICY_GRADIENT_H
#define POLICY_GRADIENT_H

#include "PolicyEvaluation.h"
#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "real.h"
#include <vector>

/** Analyitical policy gradient algorithm.
	
	This algorithm works only for discrete MDPs
 */
class PolicyGradient
{
public:
    PolicyEvaluation evaluation; ///< policy evaluation
    const DiscreteMDP* mdp; ///< the underlying MDP
    Vector starting; ///< starting state distribution
    FixedDiscretePolicy* policy; ///< policy
    std::vector<int> a_max;
    real gamma;
    int n_states;
    int n_actions;
    real Delta;
    real baseline;
    real step_size;
    PolicyGradient(const DiscreteMDP* mdp_,
                   real gamma_,
                   real step_size_);
    ~PolicyGradient();
    void Reset();
    void ModelBasedGradient(real threshold, int max_iter=-1);
    inline real getValue (int state, int action)
    {
        return evaluation.getValue(state, action);
    }
    inline real getValue (int state)
    {
        return evaluation.getValue(state);
    }
    FixedDiscretePolicy* getPolicy()
    {
        return policy;
    }
};

#endif

