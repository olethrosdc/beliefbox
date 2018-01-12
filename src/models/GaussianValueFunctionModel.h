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
#include "real.h"
#include <vector>

class GaussianValueFunctionModel : public ValueFunctionModel<Vector, int>
{
protected:
	std::vector<BayesianMultivariateRegression> model;
	int n_states;
	int n_actions;
public:
    /// Default constructor
    GaussianValueFunctionModel(int n_states_, int n_actions_) :
		n_states(n_states_),
		n_actions(n_actions_)
    {
		printf("Creating GVFM with %d states dimension, %d discrete actions\n",
			   n_states,
			   n_actions);
		for (int i=0; i<n_actions; i++) {
			model.push_back(BayesianMultivariateRegression(n_states, 1));
			model.at(i).Sampling(false);
		}
    }
    /// Default virtual destructor
    virtual ~GaussianValueFunctionModel()
    {
    }
    /// Reset the model
    virtual void Reset()
	{
		for (int i=0; i<n_actions; i++) {
			model.at(i).Reset();
		}
	}
	/// Observe a return
	virtual void AddReturnSample(const Vector& state, const int& action, const real U)
	{
		Vector Uv(1);
		Uv(0) = U;
		model.at(action).AddElement(Uv, state);
	}
	/// Calculate the values, i.e. select a model
    virtual void CalculateValues()
	{
		for (int i=0; i<n_actions; i++) {
			model.at(i).Select();
		}
	}
    /// Get the value of a state
    virtual real getValue(const Vector& state) const
	{
 		real Qmax = -INF;
		for (int i=0; i<n_actions; ++i) {
			real Qi = getValue(state, i);
			if (Qi > Qmax) {
				Qmax = Qi;
			}
		}
		return Qmax;
	}
    /// Get the value of a state-action pair
    virtual real getValue(const Vector& state, const int& action)  const
	{
		return model.at(action).generate(state)(0);
	}

};

#endif
