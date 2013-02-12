/* -*- Mode: C++; -*- */
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_NUMBER_FILE_H

#include <string>
#include <vector>
#include "RandomNumberGenerator.h"

class RandomNumberFile : RandomNumberGenerator {
protected:
    std::vector<ulong> pool;
    uint position;
public:
    RandomNumberFile(std::string filename);
    virtual ~RandomNumberFile() {}

    /// Initializes the random number generator with the computer clock.
    virtual void seed() {}
    /// Initializes the random number generator with the given long "the_seed_".
    virtual void manualSeed(ulong the_seed_)
    {
        position = the_seed_ % pool.size();
    }

    /// Returns the starting seed used.
    virtual ulong getInitialSeed()
    {
        return 0;
    }

    /// Generates a uniform 32 bits integer.
    virtual ulong random();

    /// Generates a uniform random number on [0,1[.
    virtual real uniform();
    int pool_size() {
        return pool.size();
    }
};



#endif
