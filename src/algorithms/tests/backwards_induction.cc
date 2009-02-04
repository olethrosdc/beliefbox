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

#ifdef MAKE_MAIN

#include "BackwardsInduction.h"


int main (int argc, char** argv)
{
    int n_nodes = 10;
    Vector V;
    SparseGraph G(n_nodes, true);
    MaxInduction J;
    BackwardsInduction<Vector, SparseGraph, MaxInduction> induction(V, G, J);
    return 0;
}


#endif
