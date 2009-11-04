/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef EXPECTATION_MAXIMISATION_H
#define EXPECTATION_MAXIMISATION_H

#include <vector>

/**
   \ingroup StatisticsGroup
 */
/*@{*/

template <typename Model, typename X> 
class ExpectationMaximisation
{
public:
    Model& model;
    std::vector<X> observations;
    ExpectationMaximisation(Model& model_) : model(model_)
    {
        model.Reset();
    }
    
    void Observe(X x)
    {
        observations.push_back(x);
    }
    real Iterate(int n_iterations) {
        return model.ExpectationMaximisation(observations, n_iterations);
    }
    void Reset()
    {
        observations.resize(0);
    }

};

/*@}*/

#endif
