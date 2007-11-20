/* -*- Mode: C++; -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_VARIABLE_H
#define DISCRETE_VARIABLE_H

#include <vector>

class DiscreteVariable
{
 public:
    int n_values; // number of values the discrete variable can take
    DiscreteVariable(int n) {n_values = n;}
};


// maybe this should be a class rather than a typedef
typedef std::vector<DiscreteVariable> DiscreteVector;


#endif
