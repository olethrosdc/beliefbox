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

#ifndef KNN_VALUE_FUNCTION_MODEL_H
#define KNN_VALUE_FUNCTION_MODEL_H

#include "ValueFunctionModel.h"
#include "KNNRegression.h"
#include "real.h"
#include <vector>

class KNNValueFunctionModel : public ValueFunctionModel<Vector, int>
{
protected:
	std::vector<KNNRegression*> model;
	int n_states;
	int n_actions;
	int n_neighbours;
public:
    /// Default constructor
    KNNValueFunctionModel(int n_states_, int n_actions_, int n_neighbours_) :
		n_states(n_states_),
		n_actions(n_actions_),
		n_neighbours(n_neighbours_)
    {
		printf("Creating GVFM with %d states dimension, %d discrete actions\n",
			   n_states,
			   n_actions);
		for (int i=0; i<n_actions; i++) {
			model.push_back(new KNNRegression(n_states, 1));
		}
    }
    /// Default virtual destructor
	virtual ~KNNValueFunctionModel()
    {
		for (int i=0; i<n_actions; i++) {
			delete model.at(i);
		}
    }
	/// The name of the model
	virtual const char* Name()
	{
		return "kNNVFM";
	}

    /// Reset the model
    virtual void Reset()
	{
		for (int i=0; i<n_actions; i++) {
			delete model.at(i);
			model.at(i) = new KNNRegression(n_states, 1);
		}
	}
	/// Observe a return
	virtual void AddReturnSample(const Vector& state, const int& action, const real U)
	{
		Vector Uv(1);
		Uv(0) = U;
		model.at(action)->AddElement(PointPair(state, Uv));
	}
	/// Calculate the values, i.e. select a model
    virtual void CalculateValues()
	{
	}
    /// Get the value of a state
    virtual real getValue(const Vector& state) const
	{
 		real Qmax = -INF;
		for (int i=0; i<n_actions; ++i) {
			real Qi = getValue(state, i);
			//printf ("%f ", Qi);
			if (Qi > Qmax) {
				Qmax = Qi;
			}
		}
		return Qmax;
	}
    /// Get the value of a state-action pair
    virtual real getValue(const Vector& state, const int& action)  const
	{
		Vector Uv(1);
		model.at(action)->Evaluate(state, Uv, n_neighbours);
		return Uv(0);
	}

};

#endif
