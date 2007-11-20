/* -*- Mode: C++; -*- */
/* $Revision$ */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef VECTOR_SET_H
#define VECTOR_SET_H

#include "Vector.h"
#include <list>



class VectorSet
{
 protected:
    std::list<Vector> points;
    int n;
 public:
    VectorSet(int n);
    virtual ~VectorSet();
    virtual void Add(Vector x);

};

class ConvexHull
{
protected:
    VectorSet vector_set;
public:
    virtual ~ConvexHull();
    virtual bool Add(Vector x);
    virtual bool IsInConvexHull(Vector x);
};


#endif
