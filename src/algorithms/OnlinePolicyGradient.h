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

#ifndef ONLINE_POLICY_GRADIENT_H
#define ONLINE_POLICY_GRADIENT_H

#include "PolicyEvaluation.h"
#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ContinuousPolicy.h"
#include "BasisSet.h"
#include "real.h"
#include "Sarsa.h"
#include "FeatureTD.h"
#include <vector>

/** Online policy gradient algorithms using an actor-critic architecture.
	
	Here simple cartesian features are used, together with a Sarsa implementation for a critic.

 */
class PolicyGradientActorCritic : public OnlineAlgorithm<int, int>
{
protected:
	int n_states;
	int n_actions;
	real gamma;
	real step_size;
	Sarsa critic;
	FixedDiscretePolicy policy;
	Matrix params;
	Matrix Q;
	int state, action; ///< last state and action
public:
	PolicyGradientActorCritic(int n_states_, int n_actions_, real gamma_=0.95, real step_size_ = 0.1);
	virtual ~PolicyGradientActorCritic()
	{
	}
    virtual void Reset()
	{
		state = -1;
		action = -1;
		critic.Reset();
		policy.Reset();
	}
    /// Update the actor and critic
    virtual real Observe (real reward, const int& next_state, const int& next_action);
    /// Get an action using the current policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, const int& next_state);

    virtual real getValue (const int& state, const int& action)
    {
        return critic.getValue(state, action);
    }
	virtual real getValue (const int& state)
    {
        return critic.getValue(state);
    }
	/// Update the policy for a given state and action/
	///
	/// s, a: state-action pair observed
	/// returns the gradient norm
	real GradientUpdate(int s, int a);

	void UpdatePolicy();
};

/** PGAC with features, for the continuous case.

 */
class PolicyGradientActorCriticPhi : public OnlineAlgorithm<int, Vector>
{
protected:
	BasisSet<Vector, int>& basis;
	int n_states;
	int n_actions;
	real gamma;
	real step_size;
	FeatureTD critic;
	FixedContinuousPolicy policy;
	Vector params;
	Vector state;
	int action; ///< last state and action
	bool valid_state;
public:
	PolicyGradientActorCriticPhi(BasisSet<Vector, int>& basis_,
								 int n_states_,
								 int n_actions_,
								 real gamma_=0.95,
								 real step_size_ = 0.1);
	virtual ~PolicyGradientActorCriticPhi()
	{
	}
    virtual void Reset()
	{
		valid_state = false;
		action = -1;
		critic.Reset();
		policy.Reset();
	}
    /// Update the actor and critic
    virtual real Observe (real reward, const Vector& next_state, const int& next_action);
    /// Get an action using the current policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, const Vector& next_state);

    virtual real getValue (const Vector& state, const int&  action) const
    {
        return critic.getValue(state, action);
    }
	virtual real getValue (const Vector& state) const
    {
        return critic.getValue(state);
    }
	/// Update the policy for a given state and action/
	///
	/// s, a: state-action pair observed
	/// returns the gradient norm
	real GradientUpdate(const Vector& s, int a);

	void UpdatePolicy();

};


#endif

