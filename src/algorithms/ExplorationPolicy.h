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

#ifndef EXPLORATION_POLICY_H
#define EXPLORATION_POLICY_H

#include <cassert>
#include "Matrix.h"

class ExplorationPolicy
{
public:
    virtual ~ExplorationPolicy()
    {}
    virtual int getAction(Vector Q& v) = 0;
    virtual int getAction(Matrix Q& v, int state) = 0;
};

class EpsilonGreedy
{
protected:
    real epsilon;
public:
    ExplorationPolicy(real epsilon_) :
        epsilon(epsilon_)
    {
        assert(epsilon >= 0 && epsilon <= 1);
    }
    real getEpsilon()
    {
        return epsilon;
    }
    real setEpsilon(real epsilon_)
    {
        epsilon = epsilon_;
        assert(epsilon >= 0 && epsilon <= 1);
    }
    virtual int getAction(Matrix Q& m) 
    {
    }
    virtual int getAction(Matrix Q& v, int state) 
    {
    }
    virtual ~ExplorationPolicy()
    {
    }
};
#endif
