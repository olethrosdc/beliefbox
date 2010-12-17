/* -*- Mode: c++ -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MatrixNorm.h"
#include "MathFunctions.h"

real MaxNorm(const Matrix& X)
{
    real max = 0;
    for (int i=0; i<X.Rows(); ++i) {
        for (int j=0; j<X.Columns(); ++j) {
            real f = fabs(X(i,j));
            if (f > max) {
                max = f;
            }
        }
    }
    return max;
}

real FrobeniusNorm(const Matrix& X)
{
    real norm = 0;
    for (int i=0; i<X.Rows(); ++i) {
        for (int j=0; j<X.Columns(); ++j) {
            real f = fabs(X(i,j));
            norm += f*f;
        }
    }
    return sqrt(norm);
}

real PNorm(const Matrix& X, const real p)
{
    real log_norm = LOG_ZERO;
    const real log_p = log(p);
    for (int i=0; i<X.Rows(); ++i) {
        for (int j=0; j<X.Columns(); ++j) {
            real f = fabs(X(i,j));
            log_norm = logAdd(log_norm, log_p * log(f));
        }
    }
    return exp(log_norm / p);
}
