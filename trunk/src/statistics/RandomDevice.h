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

/** This is a 'high quality' random device.

    It differs from RandomSourceRNG in that the file used is always
    opened and closed. This may be detrimental to performance.
 */
class RandomDevice : public RandomNumberGenerator 
{
protected:
    bool blocking;
public:
    RandomDevice(bool blocking_) : blocking(blocking_)
    {
    }
    virtual ~RandomDevice() 
    {
    }

    /// Nothing here
    virtual void seed()
    {
    }

    /// Nothing to be done here
    virtual void manualSeed(unsigned long the_seed_)
    {
    }

    /// Returns the starting seed used.
    virtual unsigned long getInitialSeed()
    {
        return 0;
    }
	
    /// Generates a uniform 32 bits integer.
    virtual unsigned long random()
    {
        return true_random_bits(blocking);
    }

    /// Generates a uniform random number on [0,1[.
    virtual real uniform()
    {
        return true_random(blocking);
    }
};

#endif



