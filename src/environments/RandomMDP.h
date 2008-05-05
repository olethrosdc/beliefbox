// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RANDOM_MDP_H
#define RANDOM_MDP_H

#include "DiscreteMDP.h"
#include <string>
#include <vector>



class RandomMDP {
public:
    RandomMDP(int n_actions,
              int n_states,
              real randomness);
    const DiscreteMDP* getMDP()
    {
        return mdp;
    }
 protected:
    DiscreteMDP* mdp;
    real** transitions;
    real* P_data;
    Distribution** rewards;
};

#endif
