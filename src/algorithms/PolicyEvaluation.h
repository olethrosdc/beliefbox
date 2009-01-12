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

#ifndef POLICY_EVALUATION_H
#define POLICY_EVALUATION_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "real.h"
#include <vector>

class PolicyEvaluation
{
public:
    DiscretePolicy* policy;
    const DiscreteMDP* mdp;
    real gamma;
    int n_states;
    int n_actions;
    std::vector<real> V;
    real Delta;
    real baseline;
    PolicyEvaluation(DiscretePolicy* policy_,
                     const DiscreteMDP* mdp_,
                     real gamma_,
                     real baseline_ = 0.0);
    virtual ~PolicyEvaluation();
    virtual void ComputeStateValues(real threshold, int max_iter=-1);
    inline void SetPolicy(DiscretePolicy* policy_)
    {
        policy = policy_;
    }
    void Reset();
    real getValue (int state, int action) const;
    inline real getValue (int state) const
    {
        return V[state];
    }
};

#endif

