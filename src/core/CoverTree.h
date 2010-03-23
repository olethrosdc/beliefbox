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

#ifndef COVER_TREE_H
#define COVER_TREE_H

#include "real.h"

/** A cover tree.

    Implements the algorithm "Cover trees for Nearest Neighbor",
    Beygelzimer, Kakade, Langford. Based on the technical report "Fast
    Nearest Neighbors" of Thomas Kollar.
    
 */
template <typename X>
class CoverTree
{
    struct Node
    {
        int level;
        real log_distance;
        real distance;
        X point;
        std::list<Node*> children;
        Node (int level, X& point_, real log_c)
            : level(level_),
              point(point_)
        {
            log_distance = level * log_c;
            distance = exp(log_distance);
        }
        bool Contains(X& query)
        {
            if (metric(query, point) <= distance) {
                return true;
            }
            return false;
        }
    };

    struct CoverSet
    {
        std::set<X>
        int level;
        real distance;
    };

    bool Insert(X p, CoverSet Q_i, int level)
    {
        
    }
    
};


#endif
