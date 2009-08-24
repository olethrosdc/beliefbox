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
    DiscreteVariable() {n_values = 2;}
};


// maybe this should be a class rather than a typedef
class DiscreteVector
{
public:
    std::vector<DiscreteVariable> V;
    DiscreteVector(std::vector<int>& x) : V(x.size())
    {
        for (uint i=0; i<V.size(); ++i) {
            V[i].n_values = x[i];
        }
    }
    int size(int i)
    {
        return V[i].n_values;
    }
    int size()
    {
        return V.size();
    }
    bool permute(std::vector<int>&x)
    {
        bool carry = false;
        int index = 0;
        do {
            carry = false;
            if (++x[index] == size(index)) {
                x[index] = 0;
                carry = true;
                index++;
            }
        } while (carry && index < size());

        return carry;
    }
};



#endif
