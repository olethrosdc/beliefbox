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

#ifndef KERNEL_CONDITIONAL_DENSITY_ESTIMATOR_H
#define KERNEL_CONDITIONAL_DENSITY_ESTIMATOR_H


#include "KernelDensityEstimator.h"

/** Kernel method for conditional density estimation.

    Vectors are of the form \f$(x, y)\f$.

    The conditional estimator estimates:
    \f[
    P(y | x) = P(y, x) / P(x),
    \f]
    via separate estimation of \f$P(y, x)\f$ and \f$P(x)\f$.
    
    This is a simple wrapper over KernelDensityEstimator.
 */
class KernelConditionalDensityEstimator
{
protected:
    int n_x;
    int n_y;
    KernelDensityEstimator p_xy;
    KernelDensityEstimator p_x;
public:
    KernelConditionalDensityEstimator(int n_x_dimensions,
                                      int n_y_dimensions,
                                      real initial_bandwidth,
                                      int knn);
    Vector Concatenate(const Vector& x, const Vector& y)
    {
        Vector z(x.Size() + y.Size());
        for (int i=0; i<x.Size(); ++i) {
            z(i) = x(i);
        }
        for (int i=0; i<y.Size(); ++i) {
            z(i + n_x) = y(i);
        }
        return z;
    }
    real Observe(const Vector& x, const Vector& y); 
    real pdf(const Vector& x, const Vector& y);
    real log_pdf(const Vector& x, const Vector& y);
    void BootstrapBandwidth()
    {
        p_xy.BootstrapBandwidth();
        p_x.b = p_xy.b;
        //p_x.BootstrapBandwidth();

    }
    void Show()
    {
		printf ("# Kernel CDE\n");
    }
};

#endif
