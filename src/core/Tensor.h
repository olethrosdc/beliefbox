/* -*- Mode: c++ -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TENSOR_H
#define TENSOR_H

#include <vector>
#include <cstdlib>
#include <cassert>

#include "real.h"

/** A simple tensor class.

   Access to the tensor is through an std-vector array, to obtain a
   value you need to use:

   tensor.Y(x),

   where x is std::vector<int>.

   This returns a single real value.
*/   
class Tensor
{
protected:
    std::vector<int> dimensions; ///< vector of dimensions
    real* data; ///< data storage (currently naive)
    int n; ///< number of indices
    long K; ///< total number of elements that can be stored
public:
	/// Create a tensor with the specified number of dimensions
    Tensor(std::vector<int>& dimensions_) : dimensions(dimensions_), n(dimensions.size())
    {
        K = 1;
        for (int i=0; i<n; ++i) {
            if (dimensions[i] == 0) {
                fprintf(stderr, "Tensor.h: Cannot allocate a dimension (%d) of size zero\n", i);
                exit(-1);
            }
            K *= dimensions[i];
        }
        //printf ("Making Tensor with %ld elements\n", K);
        data = (real*) calloc(K, sizeof(real));
    }

    ~Tensor()
    {
        free(data);
    }

    /// Get the value at x
    real& Y(std::vector<int>& x)
    {
        assert((int) x.size() == n);
        long index = 0;
        long M = 1;
        for (int i=0; i<n; ++i) {
            assert(x[i] >= 0 && x[i] < dimensions[i]);
            index += (long) x[i] * M;
            M *= dimensions[i];
        }
        assert(index >= 0 && index < K);
        return data[index];
    }

	/// Clear the data
    void Clear()
    {
        for (long i=0; i<K; ++i) {
            data[i] = 0;
        }
    }
};

#endif
