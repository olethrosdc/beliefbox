// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "BPSRModel.h"

int main(void)
{
    std::vector<real> rewards(4);
    rewards[0] = -1.0;
    rewards[1] = -0.1;
    rewards[2] = 0.0;
    rewards[3] = 1.0;
    int n_actions = 4;
    int n_observations = 256;
    
}


#endif
