/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_DEVICE_H
#define RANDOM_DEVICE_H

#include "real.h"
#include "RandomNumberGenerator.h"

class RandomDevice : public RandomNumberGenerator  {
protected:
    bool blocking;
public:
    virtual ~RandomDevice() {};
    /// Initializes the random number generator with the computer clock.
    static void seed();
    /// Initializes the random number generator with the given long "the_seed_".
    static void manualSeed(unsigned long the_seed_);

    /// Returns the starting seed used.
    static unsigned long getInitialSeed();
	
    /// Generates a uniform 32 bits integer.
    static unsigned long random();

    /// Generates a uniform random number on [0,1[.
    static real uniform();
};

#endif
