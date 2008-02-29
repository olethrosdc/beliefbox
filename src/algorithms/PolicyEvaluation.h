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
    DiscreteMDP* mdp;
    DiscretePolicy* p;
    real gamma;
    int n_states;
    int n_actions;
    std::vector<real> V;
    std::vector<real*> Q;
    std::vector<real> Q_data;
    real Delta;
    real baseline;
    PolicyEvaluation(DiscretePolicy* p,
		     DiscreteMDP* mdp,
		     real gamma,
		     real baseline);
    ~PolicyEvaluation();
    void Reset();
    void ComputeStateValues(real threshold, int max_iter=-1);
    void ComputeStateActionValues(real threshold, int max_iter=-1);
    real getValue (int state, int action);
    real getValue (int state);
};

#endif

