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

/** Kernel method for density estimation.

    Vectors are of the form \f$x = (y, z)\f$.

    The conditional estimator estimates directly:
    \f[
    P(y | z) = P(y, z) / P(z).
    \f]
 */
class KernelDensityEstimator
{
public:
    /// Points can be weighted.
    struct WeightedPoint
    {
        Vector x;
        real w;
        WeightedPoint(const Vector& x_, const real w_)
            : x(x_), w(w_)
        { }
    };
    KernelDensityEstimator(int n_dimensions,
                           real initial_bandwidth,
                           int knn);
    int n; ///< The number of dimensions
    real b; ///< The bandwidth
    bool change_b; ///< Whether be should be able to change
    std::list<WeightedPoint> points; ///< A list of weighted points
    real Observe(const Vector& x); 
    void AddPoint(const Vector& x, real w = 1);
    /// Return the pdf at x
    real pdf(const Vector& x)
    {
        return exp(log_pdf(x));
    }
    real log_pdf(const Vector& x);
    void BootstrapBandwidth();
    void Show()
    {
    }
protected:
    int nearest_neighbour_size; ///< what size to use for the nearest neighbour
    KDTree<WeightedPoint> kd_tree; ///< The tree, for faster access
    
    
};



#endif
