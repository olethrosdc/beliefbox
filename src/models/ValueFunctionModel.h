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

#ifndef VALUE_FUNCTION_MODEL_H
#define VALUE_FUNCTION_MODEL_H

template <typename S, typename A>
class ValueFunctionModel
{
public:
    /// Default constructor
    ValueFunctionModel()
    {
    }
    /// Default virtual destructor
    virtual ~ValueFunctionModel()
    {
    }
    /// Reset the model
    virtual void Reset() = 0;
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
