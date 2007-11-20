/* -*- Mode: C++; -*- */
// VER: $Id: Distribution.c,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Random.h"
#include "SmartAssert.h"
#include "ranlib.h"
#include "MersenneTwister.h"


void setRandomSeed(unsigned int seed)
{
    srand(seed);
	MersenneTwister::manualSeed(seed);
}

real urandom2()
{
    return MersenneTwister::uniform();
}

real urandom()
{
    real x;
    do {
        x = MersenneTwister::uniform();
    } while (x>=1.0);
    return x;
}


real true_random(bool blocking)
{
	real x;
	FILE* rand_device;
	static bool warned = false;
	if (blocking) {
		rand_device = fopen ("/dev/random", "r");
	} else {
		rand_device = fopen ("/dev/urandom", "r");
	}
	if (rand_device) {
		do {
			unsigned int i;
			fread(&i, sizeof(unsigned int), 1, rand_device);
			x = ((double) i / (double) INT_MAX);
		} while (x>=1.0);
		fclose (rand_device);
		return x;
	} else if (!warned) {
		fprintf (stderr, "Warning: true random device not available.  Using random() instead\n");
		warned = true;
	}
	
	return urandom();
}



real urandom(real min, real max)
{
    return min + ((max-min)*urandom());
}
