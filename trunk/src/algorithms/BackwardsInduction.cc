// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "BackwardsInduction.h"

#include <list>

BackwardsInduction::BackwardsInduction(Vector& N_, Graph& G_)
    : N(N_), G(G_)
{}

void BackwardsInduction::calculate()
{
    if G.hasCycles() {
        throw std::domain_error("Graph has loops, cannot use run\n");
    }
    std::vector<list> T = N;
    
    while(1) {
        
        
    }
}
