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
#include "real.h"

/** A discrete variable
    
    This variable can take \f$n\f$ different values.
 */
class DiscreteVariable
{
 public:
    int n_values; ///< number of values the discrete variable can take
    DiscreteVariable(int n) {n_values = n;}
    DiscreteVariable() {n_values = 2;}
};


/** A discrete vector
    
    A discrete vector is made up of \f$K\f$ discrete variables.
    The i-th variable takes \f$n_i\f$ values.
 */
class DiscreteVector
{
protected:
    std::vector<DiscreteVariable> V;
    int n_combinations; ///< Number of combinations
public:
        /// Construct a discrete vector from a vector of sizes
    DiscreteVector(std::vector<int>& x) : V(x.size())
    {
        n_combinations = 1;
        for (uint i=0; i<V.size(); ++i) {
            V[i].n_values = x[i];
            n_combinations *= x[i];
        }
    }
        /// Get number of combinations
    int getNCombinations() const
    {
        return n_combinations;
    }
        /// Get the number of values of the i-th variable
    int size(int i) const
    {
        return V[i].n_values;
    }
        /// Get the number of variables
    int size() const
    {
        return V.size();
    }
        /// Increase the value of one variable, carrying over to the next.
        /// Return true when the process loops.
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
    
    int getIndex(std::vector<int>& x)
    {
        int index = 0;
        int F = 1;
        for (int k=0; k<size(); ++k) {
            index += F * x[k];
            F *= size(k);
        }
        return index;
    }

};



#endif
