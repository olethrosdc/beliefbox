/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DOUBLE_KERNEL_CDE_H
#define DOUBLE_KERNEL_CDE_H


#include <vector>
#include "Vector.h"

/** Double Kernel method for conditional density estimation.

    Estimate P(x | y) = sum_c K_c(x) P(c | y).

 */
class DoubleKernelCDE
{
protected:
    int n_x; ///< dimensionality of X
    int n_y; ///< dimensionality of Y
    struct PointPair
    {
        Vector x;
        Vector y;
        PointPair(const Vector& x_,
                  const Vector& y_)
            : x(x_), y(y_)
        { 
        }
    };
    std::vector<PointPair> D;
    real b_x;
    real b_y;
public:

    DoubleKernelCDE(int n_x_dimensions,
                    int n_y_dimensions,
                    real initial_bandwidth);
    void AddPoint(const Vector& x, const Vector& y)
    {
        D.push_back(PointPair(x, y));
    }
    real Observe(const Vector& x, const Vector& y); 
    real pdf(const Vector& x, const Vector& y);
    real log_pdf(const Vector& x, const Vector& y);
    void BootstrapBandwidth(bool stochastic = false);
    void Show()
    {
		printf ("# Double Kernel CDE\n");
    }
};

#endif
