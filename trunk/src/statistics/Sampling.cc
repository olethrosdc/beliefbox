/* -*- Mode: C++; -*- */
/* VER: $Id: Sampling.c,v 1.2 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "Sampling.h"
#include <cassert>

int PropSample (std::vector<real>& w)
{
    int n = w.size();
    assert(n > 0);
    real X = UniformSample();
    real s = 0.0;
    for (int i=0; i<n; i++) {
        s += w[i];
        if (X < s) {
            return i;
        }
    }
    return rand()%n;
}
