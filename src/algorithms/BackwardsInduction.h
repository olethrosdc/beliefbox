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

#ifndef BACKWARDS_INDUCTION_H
#define BACKWARDS_INDUCTION_H


#include "Graph.h"
#include "Vector.h"

/** The backwards induction algorithm.  It needs a sort of Graph.
    This precludes it from working in anything other than a
    discrete space.
*/
//template <typename StateType, typename ActionType>
class BackwardsInduction
{
protected:
    Graph& G;
    std::list<int>& N;
public:
    // 
    BackwardsInduction(std::list<int>& N_, Graph& G_);
    void calculate(); ///< propagate values starting from terminal nodes
    void calculate_loopy(); ///< propagate values starting from any node
};

#endif

