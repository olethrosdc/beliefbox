// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "PolicyBelief.h"

real DirichletProductPolicyBelief::Update(int state, int action)
{
    assert(state >= 0 && state < n_states);
    assert(action >= 0 && action < n_actions);

    return P[state].Observe(action);
}

