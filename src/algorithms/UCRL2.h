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

#ifndef UCRL2_H
#define UCRL2_H

#include "OptimisticValueIteration.h"
#include "Matrix.h"
#include "real.h"
#include <vector>

class UCRL2 : public OnlineAlgorithms
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    int state; ///< current state

public:
    UCRL2(int n_states_,
          int n_actions_);
    ~UCRL2();
    void Reset();
    real Observe (int action, int next_state, real reward);
    real getValue (int s, int a);
};

#endif

