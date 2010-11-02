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

#ifndef ORDERED_FIXED_LIST_H
#define ORDERED_FIXED_LIST_H

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

    OrderedFixedList(const uint N_) : N(N_)
    {
        lower_bound = -REAL_RANGE;
        upper_bound = REAL_RANGE;
    }
    /** Add a (value, object) pair to the list if possible.
        
        The object is only added if either:
        
        a) There is enough space in the list.
        b) The value is smaller than the currently largest value.
    */
    const bool AddPerhaps(real x, T* object)
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
    const real UpperBound() const
    {
        return upper_bound;
    }
    const real LowerBound() const
    {
        return lower_bound;
    }
    const int max_size() const
    {
        return S.size();
    }
    const int size() const
    {
        return N;
    }

};


#endif
