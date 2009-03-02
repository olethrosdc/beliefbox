/* -*- Mode: C++; -*- */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

#include "real.h"
#include <cmath>

/// Maybe this class is a bit useless.
class  RandomNumberGenerator {
public:
    virtual ~RandomNumberGenerator() {}

    /// Initializes the random number generator with the computer clock.
    virtual void seed() = 0;
    /// Initializes the random number generator with the given long "the_seed_".
    virtual void manualSeed(unsigned long the_seed_) = 0;

    /// Returns the starting seed used.
    virtual unsigned long getInitialSeed()  = 0;
	
    /// Generates a uniform 32 bits integer.
    virtual unsigned long random() = 0;

    /// Generates a uniform random number in [0,1[.
    virtual real uniform() = 0;

    /// Generates a uniform random number in [0,n)
    virtual int discrete_uniform(int n)
    {
        return (int) floor(uniform()*((real) n));
    }
};

#endif
