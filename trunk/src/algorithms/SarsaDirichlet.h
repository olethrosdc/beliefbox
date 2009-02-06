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

#ifndef SARSA_DIRICHLET_H
#define SARSA_DIRICHLET_H

#include "Sarsa.h"

class SarsaDirichlet : public Sarsa
{
protected:
    Matrix C; ///< counts
public:
    SarsaDirichlet(int n_states_,
          int n_actions_,
          real gamma_,
          real lambda_,
          real alpha_,
          VFExplorationPolicy* exploration_policy_,
          real initial_value_= 0.0,
          real baseline_ = 0.0);
    virtual ~SarsaDirichlet();
    /// Partial SARSA observation (can be used with eligibility traces)
    virtual real Observe (real reward, int next_state, int next_action);
    /// Get an action using the current exploration policy.
    virtual int Act(real reward, int next_state);
};

#endif

