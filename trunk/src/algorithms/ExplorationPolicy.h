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
#include "Random.h"

class ExplorationPolicy
{
public:
    virtual ~ExplorationPolicy()
    {}
    virtual int getAction(Vector& Q) = 0;
    virtual int getAction(Matrix& Q, int state) = 0;
};

class EpsilonGreedy : public ExplorationPolicy
{
protected:
    int n_actions;
    real epsilon;
public:
    EpsilonGreedy(int n_actions_, real epsilon_) :
        n_actions(n_actions_), epsilon(epsilon_)
    {
	assert(n_actions > 0);
        assert(epsilon >= 0 && epsilon <= 1);
    }

    virtual ~EpsilonGreedy()
    {
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
    virtual int getAction(Vector& Q) 
    {
	if (urandom() < epsilon) {
	    return (int) floor(urandom(0, n_actions));
	}
	return ArgMax(&Q);
    }

    virtual int getAction(Matrix& Q, int state) 
    {
	if (urandom() < epsilon) {
	    return (int) floor(urandom(0.0, n_actions));
	}
	int argmax = 0;
	real max = Q(state, argmax);
	for (int a=1; a<n_actions; ++a) {
	    real Qsa = Q(state, a);
	    if (Qsa > max) {
		max = Qsa;
		argmax = a;
	    }
	}
	return argmax;
    }
};
#endif
