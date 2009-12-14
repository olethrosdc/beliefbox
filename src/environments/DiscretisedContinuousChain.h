// -*- Mode: c++ -*-
// copyright (c) 2008-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// copyright (c) 2003-2008 Michail G. Lagoudakis
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETISED_CONTINUOUS_CHAIN_H
#define DISCRETISED_CONTINUOUS_CHAIN_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

class DiscretisedContinuousChain : public Environment<int, int>
{
protected:
    void Simulate();
    bool endsim;
    real position;
    void DiscretiseState()
    {
        if (endsim) {
            state = 0;
            return;
        } 
        real diff = (position + 1.0) / 2.0;
        state = 1 + (int) round(diff*(n_states - 2));
        if (state <= 0) state = 1;
        else if (state >= (int) n_states) state = n_states - 1;
    }
public:
    DiscretisedContinuousChain(int n_states);
    virtual ~DiscretisedContinuousChain();
    virtual void Reset();
    virtual bool Act(int action);
    virtual void Simulate(int action);
};








#endif
