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

#ifndef KNNMODEL_H
#define KNNMODEL_H


template<typename X>
class KNN
{
protected:
    std::vector<X> point_set;
public:	
    struct PointDistance
    {
        X& x;
        real d;
    };

    void AddElement(X x)
    {
        point_set.push_back(x);
    }

    std::vector<PointDistance> FindKNearestNeigbours(X x)
    {
        std::vector<real> d(point_set.size());
        for (uint i=1; i<point_set.size(); ++i) {
            x.distance(L.)
        }
    }
};


#endif
