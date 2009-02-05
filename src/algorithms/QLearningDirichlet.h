// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef Q_LEARNING_DIRICHLET_H
#define Q_LEARNING_DIRICHLET_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "real.h"
#include "QLearning.h"
#include <vector>

class QLearningDirichlet : public QLearning
{
protected:
    Matrix C; ///< matrix of counts
public:
    QLearningDirichlet(int n_states_,
                       int n_actions_,
                       real gamma_,
                       real lambda_,
                       real alpha_,
                       VFExplorationPolicy* exploration_policy_,
                       real initial_value_= 0.0,
                       real baseline_ = 0.0);
    virtual ~QLearningDirichlet()
    {
    }
    virtual real Observe (int action, int next_state, real reward);
};

#endif

