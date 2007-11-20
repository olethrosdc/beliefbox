/* -*- Mode: C++; -*- */
/* VER: $Id: Sampling.h,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SAMPLING_H
#define SAMPLING_H

#include <cstdlib>
#include <vector>
#include "real.h"


inline real UniformSample()
{
    return drand48();
}

int PropSample (std::vector<real>& w);

#endif
