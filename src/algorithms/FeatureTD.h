// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FEATURE_TD_H
#define FEATURE_TD_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "BasisSet.h"
#include "Critic.h"
#include <vector>

class FeatureTD : public Critic<Vector, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real alpha; ///< learning rate
	BasisSet<Vector, int>& features;
    Vector params; ///< parameters
    Vector state; ///< current state
    int action; ///< current action
	bool valid_state; ///< are we at the start of an episode?
public:
    FeatureTD(int n_states_,
			  int n_actions_,
			  BasisSet<Vector, int>& features_,
			  real gamma_,
			  real alpha_=0.5);
    virtual ~FeatureTD();
    virtual void Reset();
    /// Partial observation
    virtual real Observe (real reward,
						  const Vector& next_state,
						  const int& next_action);
    virtual real getValue (const Vector& state, const int& action) const;
	virtual real getValue (const Vector& state) const;

};

#endif

