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

#ifndef BELIEF_H
#define BELIEF_H

#include "Object.h"

/** This is a generic belief container class.
	
	The purpose of this class is to serve as a container for various
	types of beliefs; both those arising from simple conjugate priors
	of i.i.d. observations, as well as more complex ones.
	
	For this reason, the implementation of a belief object should be
	relatively simple, with most of the functionality being
	implemented in the contained classes.
	
 */
class Belief : public Object
{
 public:
    virtual ~Belief() {}	
};

#endif
