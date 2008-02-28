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

FixedDiscretePolicy::FixedDiscretePolicy (Vector p)
{
    this->p = p;
    assert (fabs(this->p.Sum() - 1.0) < 0.00001);
}

FixedDiscretePolicy::~FixedDiscretePolicy()
{
}

int FixedDiscretePolicy::SelectAction()
{
    int n = p.Size();
    real x = urandom();
    real s = 0.0;
    for (int a=0; a<n; ++a) {
	s += p[a];
	if (s>x) {
	    return a;
	}
    }
    return n-1;
}

void FixedDiscretePolicy::Observe (int previous_state, int action, real r, int next_state)
{
}

void FixedDiscretePolicy::Reset()
{
}
