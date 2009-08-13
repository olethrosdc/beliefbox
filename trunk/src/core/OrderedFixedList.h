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
class OrderedFixedList
{
protected:
    const uint N;
    real lower_bound;
    real upper_bound;
public:
    std::list<real> S;
    OrderedFixedList(uint N_) : N(N_)
    {
        lower_bound = -RAND_MAX;
        upper_bound = RAND_MAX;
    }
    bool AddPerhaps(real x)
    {
        std::list<real>::iterator it;
        if (S.size() < N) {
            S.push_back(x);
            S.sort();
            lower_bound = S.front();
            upper_bound = S.back();
            return true;
        }
        if (x < upper_bound) {
            S.pop_back();
            S.push_back(x);
            S.sort();
            upper_bound = S.back();
            return true;
        }
        return false;
    }

    

};


#endif
