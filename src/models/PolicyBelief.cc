// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "PolicyBelief.h"
#include "DiscretePolicy.h"
#include "Demonstrations.h"


/// Update distribution from a single state-action observation
real DirichletProductPolicyBelief::Update(int state, int action)
{
    assert(state >= 0 && state < n_states);
    assert(action >= 0 && action < n_actions);

    return P[state].Observe(action);
}

/// Create a new policy by sampling from the current distribution
DiscretePolicy* DirichletProductPolicyBelief::Sample() const
{
    std::vector<Vector> p_sa(n_states);
    for (int i=0; i<n_states; ++i) {
        p_sa[i] = P[i].generate();
    }
    DiscretePolicy* policy = new FixedDiscretePolicy(p_sa);
    return policy;
}

/// Create a new policy, the same as the expected policy. 
DiscretePolicy* DirichletProductPolicyBelief::getExpectedPolicy() const
{
    std::vector<Vector> p_sa(n_states);
    for (int i=0; i<n_states; ++i) {
        p_sa[i] = P[i].generate();
    }
    DiscretePolicy* policy = new FixedDiscretePolicy(p_sa);
    return policy;
}

/// Get the posterior
real DirichletProductPolicyBelief::CalculatePosterior(Demonstrations<int, int>& D)
{
    real log_p = 0;
    for (uint k=0; k<D.trajectories.size(); ++k) {
        for (uint t=0; t<D.trajectories[k].x.size(); ++t) {
            int s = D.trajectories[k].x[t].first;
            int a = D.trajectories[k].x[t].second;
            assert(s >= 0 && s < n_states);
            assert(a >= 0 && a < n_actions);
            log_p += log(Update(s, a));
        }
    }
    return exp(log_p);
}


real DirichletProductPolicyBelief::getLogDensity(const DiscretePolicy& policy) const
{
	real log_pdf = 0;
	for (int s=0; s<n_states; ++s) {
		log_pdf += P[s].log_pdf(policy.getActionProbabilities(s));
	}
	return log_pdf;
}
real DirichletProductPolicyBelief::getDensity(const DiscretePolicy& policy) const
{
	return exp(getLogDensity(policy));
}


