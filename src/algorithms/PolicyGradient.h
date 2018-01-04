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

#ifndef POLICY_ITERATION_H
#define POLICY_ITERATION_H

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
    PolicyEvaluation* evaluation; ///< policy evaluation
    const DiscreteMDP* mdp;
	Vector starting; ///< starting state distribution
    FixedDiscretePolicy* policy;
    std::vector<int> a_max;
    real gamma;
    int n_states;
    int n_actions;
    real Delta;
    real baseline;
    PolicyGradient(PolicyEvaluation* evaluation_,
				   const DiscreteMDP* mdp_,
				   const Vector& starting_,
				   real gamma_);
    ~PolicyGradient();
    void Reset();
    void ComputeStateValues(real threshold, int max_iter=-1);
    inline real getValue (int state, int action)
    {
        return evaluation->getValue(state, action);
    }
    inline real getValue (int state)
    {
        return evaluation->getValue(state);
    }
};

#endif

