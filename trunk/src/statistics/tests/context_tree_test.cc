/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "ContextTree.h"
#include "Random.h"
#include "RandomNumberGenerator.h"
#include "MersenneTwister.h"
#include <ctime>

int main(void)
{
	int depth = 3;
	int n_symbols = 2;
	ContextTree tree(n_symbols, depth);
	MersenneTwisterRNG mt;
	RandomNumberGenerator* rng = &mt;
	
	int T = 10;
	for (int t=0; t<T; ++t) {
		int x = rng->discrete_uniform(n_symbols);
		tree.Observe(x, x);
	}
	return 0;
}

#endif
