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

/// Value-function-based exploration policy
///
/// Examples: epsilon-greedy, softmax
class VFExplorationPolicy
{
public:
    virtual ~VFExplorationPolicy()
    {}
    virtual void Observe(real reward, int state) = 0;
    virtual int SelectAction() = 0;
    virtual void setValueMatrix(Matrix* Q) = 0;
    virtual DiscretePolicy* getFixedPolicy() = 0;
};

class EpsilonGreedy : public VFExplorationPolicy
{
protected:
    int n_actions;
    real epsilon;
    int state;
    Matrix* Q;
public:
    EpsilonGreedy(int n_actions_, real epsilon_) :
        n_actions(n_actions_), epsilon(epsilon_), Q(NULL)
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
        return epsilon;
    }
    virtual void setValueMatrix(Matrix* Q_)
    {
        Q = Q_;
    }

    virtual int SelectAction() 
    {
	if (urandom() < epsilon) {
	    return (int) floor(urandom(0.0, n_actions));
	}
	int argmax = 0;
	real max = (*Q)(state, argmax);
	for (int a=1; a<n_actions; ++a) {
	    real Qsa = (*Q)(state, a);
	    if (Qsa > max) {
		max = Qsa;
		argmax = a;
	    }
	}
	return argmax;
    }
    virtual void Observe(real reward, int state)
    {
        this->state = state;
    }

    virtual DiscretePolicy* getFixedPolicy() 
    {
        Serror ("Not implemented\n");
        exit(-1);
        return NULL;
    }

};
#endif
