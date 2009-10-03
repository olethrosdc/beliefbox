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

class Tensor
{
protected:
    std::vector<int> X;
    real* data;
    int n;
    long K;
public:
    Tensor(std::vector<int>& X_) : X(X_), n(X.size())
    {
        K = 1;
        for (int i=0; i<n; ++i) {
            if (X[i] == 0) {
                fprintf(stderr, "Tensor.h: Cannot allocate a dimension (%d) of size zero\n", i);
                exit(-1);
            }
            K *= X[i];
        }
        //printf ("Making Tensor with %ld elements\n", K);
        data = (real*) calloc(K, sizeof(real));
    }

    ~Tensor()
    {
        free(data);
    }

    // get value
    real& Y(std::vector<int>& x)
    {
        assert((int) x.size() == n);
        long index = 0;
        long M = 1;
        for (int i=0; i<n; ++i) {
            assert(x[i] >= 0 && x[i] < X[i]);
            index += (long) x[i] * M;
            M *= X[i];
        }
        assert(index >= 0 && index < K);
        return data[index];
    }
    void Reset()
    {
        for (long i=0; i<K; ++i) {
            data[i] = 0;
        }
    }
};

#endif
