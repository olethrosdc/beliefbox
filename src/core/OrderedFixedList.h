/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ORDERED_FIXED_LIST
#define ORDERED_FIXED_LIST

#include <iostream>
#include <list>
#include "real.h"

template <typename T>
class OrderedFixedList
{
protected:
    const uint N;
    real lower_bound;
    real upper_bound;
public:
    std::list<std::pair<real, T*> > S;
    OrderedFixedList(uint N_) : N(N_)
    {
        lower_bound = -RAND_MAX;
        upper_bound = RAND_MAX;
    }
    bool AddPerhaps(real x, T* object)
    {
        std::list<real>::iterator it;
        if (S.size() < N) {
            S.push_back(std::make_pair(x, object));
            S.sort();
            lower_bound = S.front().first;
            upper_bound = S.back().first;
            return true;
        }
        if (x < upper_bound) {
            S.pop_back();
            S.push_back(std::make_pair(x, object));
            S.sort();
            upper_bound = S.back().first;
            return true;
        }
        return false;
    }
    real UpperBound()
    {
        return upper_bound;
    }
    real LowerBound()
    {
        return lower_bound;
    }

};


#endif
