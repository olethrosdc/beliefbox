/* -*- Mode: C++; -*- */
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MULTIVARIATE_NORMAL_VALUE_FUNCTION_MODEL_H
#define MULTIVARIATE_NORMAL_VALUE_FUNCTION_MODEL_H

#include "ValueFunctionModel.h"
#include "MultivariateNormalUnknownMeanPrecision.h"

template <>
class MultivariateNormalValueFunctionModel<Vector, int>
{
protected:
	std::vector<MultivariateNormalUnknownMeanPrecision*> distribution;
	int n_states;
	int n_actions;
public:
    /// Default constructor
    MultivariateNormalValueFunctionModel(int n_states_,
										 int n_actions_)
		: n_states(n_states_),
		  n_actions(n_actions_)
    {
		distribution.resize(n_actions)
		for (uint i=0; i<distribution.size(); i++) {
			distribution.at(i) = NULL;
		}
		Reset();
    }
    /// Default virtual destructor
    virtual ~MultivariateNormalValueFunctionModel()
    {
    }
	virtual const char* Name()
	{
		return "Multivariate Normal Value Function Model";
	}
    /// Reset the model
    virtual void Reset()
	{
		Vector mean(n_states);
		real tau = 1.0;
		real alpha = 1.0;
		Matrix T = Matrix::Unity(n_states, n_states);
		for (uint i=0; i<distribution.size(); i++) {
			delete distribution.at(i);
			distribution.at(i) = new MultivariateNormalUnknownMeanPrecision(mean, tau, alpha, T);
		}
	}
	/// Observe a return
	virtual void AddReturnSample(const S& state, const A& action, const real U) = 0;
	/// Calculate the values
    virtual void CalculateValues() = 0;
    /// Get the value of a state
    virtual real getValue(const S& state) const = 0;
    /// Get the value of a state-action pair
    virtual real getValue(const S& state, const A& action)  const = 0;

};

#endif
