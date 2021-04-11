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

#include "FeatureTD.h"

FeatureTD::FeatureTD(int n_states_,
					 int n_actions_,
					 BasisSet<Vector, int>& features_,
					 real gamma_,
					 real alpha_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      alpha(alpha_),
	  features(features_)
{
    assert (alpha >= 0 && alpha <= 1);
    assert (gamma >=0 && gamma <= 1);
	Reset();
}


FeatureTD::~FeatureTD() 
{
}

/// Use when starting a new episode.
void FeatureTD::Reset()
{
}

/// Observe a state-action-reward-state-action tuplet.
///
/// No eligibility trace update!
real FeatureTD::Observe (real reward,
						 const Vector& next_state,
						 const int& next_action)
{
}

int FeatureTD::Act(real reward, const Vector& next_state)
{
	return 0;
}




