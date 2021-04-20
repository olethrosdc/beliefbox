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
	  features(features_),
	  params(features.size())
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
/// We need to perform the update
/// V += a * (r - g V')

real FeatureTD::Observe (real reward,
						 const Vector& next_state,
						 const int& next_action)
{
	features.Evaluate(state);
	Vector phi_t = features.F();
	features.Evaluate(next_state);
	Vector phi_tn = features.F();
	real V = Product(phi_t, params);
	real V_n = reward + gamma * Product(phi_tn, params);
	real delta = V_n - V;
	params += alpha * phi_t * delta; // sanity check: increase when delta >0
	state = next_state;
	action = next_action;
}
real FeatureTD::getValue (const Vector& state, const int& action) const
{
	Serror("Not implemented\n");
	exit(-1);
}
real FeatureTD::getValue (const Vector& state) const
{
	features.Evaluate(state);
	Vector phi_t = features.F();
	real V = Product(phi_t, params);
	return V;
}


int FeatureTD::Act(real reward, const Vector& next_state)
{
	return 0;
}




