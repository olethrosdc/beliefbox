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
#ifndef RANDOM_SOURCE_RNG_H
#define RANDOM_SOURCE_RNG_H

#include "RandomNumberGenerator.h"
#include "Random.h"

/// This random number generator can be passed
/// as an argument.
class RandomSourceRNG : public RandomNumberGenerator
{
protected:
    FILE* rand_device;
    std::string rand_device_name;
public:
    RandomSourceRNG(bool blocking);
    virtual ~RandomSourceRNG();

    /// Initializes the random number generator with the computer clock.
    virtual void seed()
    {
    }
    /// Initializes the random number generator with the given long "the_seed_".
    virtual void manualSeed(unsigned long the_seed_)
    {
    }
    /// Returns the starting seed used.
    virtual unsigned long getInitialSeed() 
    {
        return 0;
    }
	
    /// Generates a uniform 32 bits integer.
    virtual unsigned long random();

    /// Generates a uniform random number in [0,1[.
    virtual real uniform();

};

#endif
