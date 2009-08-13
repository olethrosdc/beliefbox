// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_POMDP_H
#define DISCRETE_POMDP_H

#include "DiscreteMDP.h"

class DiscreteMDP
{
public:
    DiscreteMDP* mdp; ///< The underlying MDP
    std::vector<Distribution*> R;
    int GenererateObservation();
};

#endif
