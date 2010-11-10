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
#include "Random.h"
#include <cmath>

/// A general random number generator
class  RandomNumberGenerator {
public:
    virtual ~RandomNumberGenerator() {}

    /// Initializes the random number generator with the computer clock.
    virtual void seed() 
    {
        manualSeed(time(NULL));
    }
    /// Initializes the random number generator with the given long "the_seed_".
    virtual void manualSeed(unsigned long seed) = 0;

    /// Returns the starting seed used.
    virtual unsigned long getInitialSeed()  = 0;
	
    /// Generates a uniform 32 bits integer.
    virtual unsigned long random() = 0;

    /// Generates a uniform random number in [0,1[.
    virtual real uniform() = 0;

    /// Generates a uniform random number in [0,n)
    inline int discrete_uniform(int n)
    {
        return (int) floor(uniform()*((real) n));
    }

    inline real uniform(real lower_bound, real upper_bound)
    {
        return lower_bound + (upper_bound - lower_bound) * uniform();
    }
};

class DefaultRandomNumberGenerator : public RandomNumberGenerator
{
protected:
    ulong initial_seed;
 public:
    virtual ~DefaultRandomNumberGenerator() {}

    /// Initializes the random number generator with the given long "the_seed_".
    virtual void manualSeed(unsigned long seed) 
    {
        initial_seed = seed;
        setRandomSeed(seed);
    }

    /// Returns the starting seed used.
    virtual unsigned long getInitialSeed() 
    {
        return 0;
    }
	
    /// Generates a uniform 32 bits integer.
    virtual unsigned long random() 
    {
        return lrandom();
    }
    /// Generates a uniform random number in [0,1[.
    virtual real uniform()
    {
        return urandom();
    }
};
#endif
