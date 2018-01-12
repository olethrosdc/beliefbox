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

#ifndef GAUSSIAN_VALUE_FUNCTION_MODEL_H
#define GAUSSIAN_VALUE_FUNCTION_MODEL_H

#include "ValueFunctionModel.h"
#include "BayesianMultivariateRegression.h"
#include <vector>

template <>
class GaussianValueFunctionModel<Vector, int> : public ValueFunctionModel
{
protected:
	std::vector<BayesianMultivariateRegression> model;
	int n_states;
	int n_actions;
public:
    /// Default constructor
    ValueFunctionModel(int n_states_, int n_actions_) :
		n_states(n_states_),
		n_actions(n_actions_),
		model(n_states)
    {
		
    }
    /// Default virtual destructor
    virtual ~ValueFunctionModel()
    {
    }
    /// Reset the model
    virtual void Reset()
	{
		
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
