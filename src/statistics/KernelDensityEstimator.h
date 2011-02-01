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

#ifndef KERNEL_DENSITY_ESTIMATOR_H
#define KERNEL_DENSITY_ESTIMATOR_H

#include "NormalDistribution.h"
#include "KDTree.h"
#include "Vector.h"
#include <list>


class KernelDensityEstimator
{
public:
    struct WeightedPoint
    {
        Vector x;
        real w;
        WeightedPoint(const Vector& x_, const real w_)
            : x(x_), w(w_)
        { }
    };
    KernelDensityEstimator(int n_dimensions,
                           real initial_bandwidth);
    int n; ///< The number of dimensions
    real b; ///< The bandwidth
    bool change_b;
    std::list<WeightedPoint> points; ///< A list of weighted points
    real Observe(const Vector& x);
    void AddPoint(const Vector& x, real w = 1);
    real pdf(const Vector& x);
    real log_pdf(const Vector& x);
    void BootstrapBandwidth();
    void Show()
    {
    }
    
};



#endif
