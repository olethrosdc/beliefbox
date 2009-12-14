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
#ifndef DISCRETE_HIDDEN_MARKOV_MODEL_ONLINE_EM_H
#define DISCRETE_HIDDEN_MARKOV_MODEL_ONLINE_EM_H

#include "DiscreteHiddenMarkovModel.h"
#include "Tensor.h"
#include "Vector.h"
#include <vector>


/**
   \ingroup StatisticsGroup
 */
/*@{*/

class DiscreteHiddenMarkovModelOnlineEM  : public DiscreteHiddenMarkovModel
{
protected:
    Tensor* Phi;
    Vector q; ///< probabiltiy of being in a state at time T
    void InitPhi();
    real T;
public:
    DiscreteHiddenMarkovModelOnlineEM(Matrix& Pr_S, Matrix& Pr_X);
    DiscreteHiddenMarkovModelOnlineEM(int n_states_, int n_observations_);
    void Reset();
    virtual ~DiscreteHiddenMarkovModelOnlineEM ();
    int getCurrentState()
    {
        return current_state;
    }
    real Observe(int x);
    void updateSufficientStatistics();
    Vector getPrediction();
    int predict()
    {
        return ArgMax(getPrediction());
    }
};


/*@}*/
#endif
