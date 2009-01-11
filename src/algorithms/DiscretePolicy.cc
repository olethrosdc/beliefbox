/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscretePolicy.h"
//#include "MDPDistribution.h"

#include "Random.h"
#include <cmath>

FixedDiscretePolicy::FixedDiscretePolicy (std::vector<Vector> p)
    : DiscretePolicy()
{
    state = 0;
    this->p = p;
    for (uint i=0; i<p.size(); i++) {
        assert (fabs(this->p[i].Sum() - 1.0) < 0.00001);
    }
}

FixedDiscretePolicy::~FixedDiscretePolicy()
{
}

int FixedDiscretePolicy::SelectAction()
{
    
    int n = p[state].Size();
    real x = urandom();
    real s = 0.0;
    for (int a=0; a<n; ++a) {
        s += p[state][a];
        if (s>x) {
            return a;
        }
    }
    return n-1;
}

void FixedDiscretePolicy::Observe (int previous_state, int action, real r, int next_state)
{
    state = next_state;
}

void FixedDiscretePolicy::Observe (real r, int next_state)
{
    state = next_state;
}


void FixedDiscretePolicy::Reset()
{
    state = 0;
}

real FixedDiscretePolicy::getActionProbability(int action)
{
	return p[state][action];
}


real FixedDiscretePolicy::getActionProbability(int state, int action)
{
	return p[state][action];
}



