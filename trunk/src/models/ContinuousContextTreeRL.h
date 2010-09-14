#if 0
/* -*- Mode: c++;  -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTEXT_TREE_RL_H
#define CONTEXT_TREE_RL_H

#include <vector>
#include <list>
#include "real.h"
#include "Vector.h"
#include "Ring.h"
#include "BetaDistribution.h"


/** A continuous context tree implementation.
    
    This tree partitions an \f$nf\$-dimensional space.  It is very
    similar to ContextTreeRL, the only difference being that the
    partition is done in the state-space.

    The simplest models are:

    1. Each context model predicts only observations within that level.
    So, \f$\Pr(z_{t+1} \in c' \mid z_t \in c)\f$, with \f$c, c' \in C_k\f$.
    
    2. Each context model predicts next contexts in general
    So, \f$\Pr(z_{t+1} \in c' \mid z_t \in c)\f$, w'ith \f$c, c' \in C\f$.

    3. There is a full density estimation.
    
    @see ContextTreeRL, ConditionalKDContextTree
*/
class ContinuousContextTreeRL
{
};



#endif

#endif
