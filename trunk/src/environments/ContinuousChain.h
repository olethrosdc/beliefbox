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

#ifndef CONTINUOUS_CHAIN_H
#define CONTINUOUS_CHAIN_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"

/** A continuous chain.
    
    The state space is \f$S = [-1,1]\f$. 
    There are three actions, each increasing or decreasing the state's value.
    When the state is +1 or -1, there is a transition back.
    
*/
class ContinuousChain : public Environment<Vector, int>
{
protected:
    void Simulate();
    bool endsim;
public:
    ContinuousChain();
    virtual ~ContinuousChain();
    virtual void Reset();
    virtual bool Act(int action);
    virtual void Simulate(int action);
};








#endif
