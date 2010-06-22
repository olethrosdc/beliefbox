/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "RandomSourceRNG.h"
#include <stdexcept>

RandomSourceRNG::RandomSourceRNG(bool blocking)
{
 	if (blocking) {
		rand_device = fopen ("/dev/random", "r");
	} else {
		rand_device = fopen ("/dev/urandom", "r");
	}
    if (!rand_device) {
        throw std::runtime_error ("Could not open random device\n");
    }
}

RandomSourceRNG::~RandomSourceRNG()
{
    if (rand_device) {
        fclose(rand_device);
    }
}
/// Generates a uniform 32 bits integer.
unsigned long RandomSourceRNG::random()
{
 	unsigned long x;
    fread(&x, sizeof(unsigned long), 1, rand_device);
    return x;
}

/// Generates a uniform random number in [0,1[.
real RandomSourceRNG::uniform()
{	
    real x;
	do {
        unsigned int i;
        fread(&i, sizeof(unsigned int), 1, rand_device);
        x = ((double) i / (double) std::numeric_limits<int>::max());
    } while (x>=1.0);
    return x;
}
