/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscretePolyaTree.h"

DiscretePolyaTree::DiscretePolyaTree(DiscreteVector values_)
    : values(values_)
{
    // set up the parameters recursively
    root = new PolyaNode(values);
}

void DiscretePolyaTree::Observe(std::vector<int>& x)
{
    root->Observe(x);
}



