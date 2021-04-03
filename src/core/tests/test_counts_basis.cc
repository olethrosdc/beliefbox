/* -*- Mode: c++ -*- */
// copyright (c) 2021 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifdef MAKE_MAIN

#include "BasisSet.h"

int main(void)
{
	int n_states = 3;
	int n_actions = 2;

	CountsBasis basis(n_states, n_actions);
	basis.Reset();
	for (int i=0; i<100; ++i) {
		int s = urandom(0, n_states);
		int a = urandom(0, n_actions);
		real r = urandom();
		basis.Observe(a, r, s);
		Vector f = basis.F();
		f.print(stdout);
	}
}
#endif
