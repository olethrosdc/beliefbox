// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef Q_LEARNING_H
#define Q_LEARNING_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "real.h"
#include <vector>

class QLearning
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real lambda; ///< eligibility trace decay rate
    real alpha; ///< learning rate 
    real initial_value; ///< initial value for Q values
    real baseline; ///< baseline reward

    Matrix Q;
    Matrix el;

    int s; ///< current state

public:
    QLearning(int n_states_,
              int n_actions_,
              real gamma_,
              real lambda_,
              real alpha_ = 0.1,
              real initial_value_= 0.0,
              real baseline_ = 0.0);
    ~QLearning();
    void Reset();
    int Act (int next_state, real reward);
    real getValue (int state, int action);
    
};

#endif

