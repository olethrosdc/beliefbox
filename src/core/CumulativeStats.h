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

#ifndef CUMULATIVE_STATS_H
#define CUMULATIVE_STATS_H

#include "Vector.h"
#include "Matrix.h"

/// Cumulative statistics
struct CumulativeStats
{
    int T; ///< sequence length
    int K; ///< number of sequences
    Matrix C; ///< cumulative sum
    int current_sequence;
    CumulativeStats(int T_, int K_) : T(T_),
                                      K(K_),
                                      C(K,T),
                                      current_sequence(0)
    {
    }
    CumulativeStats(const Matrix& C_) : T(C_.Columns()),
                                        K(C_.Rows()),
                                        C(C_),
                                        current_sequence(0)
    {
    }
        
    void SetSequence(int seq)
    {
        assert (seq >= 0 && seq < K);
        current_sequence = seq;
    }
    void SetValue(int t, real x)
    {
        if (t == 0) {
            C(current_sequence, t) = x;
        } else {
            C(current_sequence, t) = x + C(current_sequence, t - 1); 
        }
    }
    void Sort()
    {
        for (int t=0; t<T; ++t) {
            C.SortColumn(t);
        }
    }
    Vector BottomPercentile(real p)
    {
        assert(p>=0 && p<=1);
        int k = (int) ceil(p * (real) K);
        k = std::min(k, C.Rows() - 1);
        assert(k >= 0 && k < K);
        return C.getRow(k);
    }
    Vector TopPercentile(real p)
    {
        assert(p>=0 && p<=1);
        int k = (int) ceil(p * (real) K);
        k = std::min(k, C.Rows() - 1) + 1;
        assert(K - k >= 0 && K - k < K);
        return C.getRow(K - k);
    }
    Vector Mean()
    {
        Vector mean(T);
        for (int k=0; k<K; ++k) {
            mean += C.getRow(k);
        }
        return mean / (real) K;
    }
};

#endif /* CUMULATIVE_STATS_H */
