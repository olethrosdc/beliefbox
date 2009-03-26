//  -*- Mode: c++ -*-
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

#ifndef DISCRETE_STATE_AGGREGATE
#define DISCRETE_STATE_AGGREGATE

#include <set>

/// A container for states
class DiscreteStateAggregate
{
public:
    std::set<int> S;
    void add(int state)
    {
        S.insert(state);
    }
    bool contains(int state) const
    {
        return (S.find(state) != S.end());
    }
    int size() const
    {
        return S.size();
    }
};

#endif
