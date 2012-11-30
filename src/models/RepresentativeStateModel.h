/// -*- Mode: c++ -*-
// copyright (c) 2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REPRESENTATIVE_STATE_MODEL_H
#define REPRESENTATIVE_STATE_MODEL_H

#include "MDP.h"

template <class Model, class S, class A>
class RepresentativeStateModel : class MDP<S, A>
{
protected:
    Model model;
    std::vector<S> states;
    Vector V;
public:
    RepresentativeStateModel(const Model& model_, const std::vector<S>& states_) :
        model(model_),
        states(states_)
    {
    }
    void AddState(const S& state)
    {
        states.push_back(state);
    }
    ComputeStateValues(real threshold, int max_iter = -1)
    {
        ValueIteration VI;
        MDP
    }
    
};

#endif
