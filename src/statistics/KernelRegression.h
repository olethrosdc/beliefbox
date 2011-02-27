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

#ifndef KERNEL_REGRESSION_H
#define KERNEL_REGRESSION_H

#include "NormalDistribution.h"
#include "KDTree.h"
#include "Vector.h"
#include <list>

/** Kernel method for regression.

 */
class KernelRegression
{
public:
    /// Points can be weighted.
    struct PointPair
    {
        Vector x;
        Vector y;
        PointPair(const Vector& x_, const Vector y_)
            : x(x_), y(y_)
        { }
    };
    KernelRegression(int n_dimensions_x,
					 int n_dimensions_y,
					 real initial_bandwidth,
					 int knn = 0);
    int n_x; ///< The number of dimensions in x
    int n_y; ///< The number of dimensions in y
    real b; ///< The bandwidth
    bool change_b; ///< Whether be should be able to change
    std::list<PointPair> points; ///< A list of weighted points
    real Observe(const Vector& x, const Vector& y); 
    void AddPoint(const Vector& x, const Vector& y);
    Vector expected_value(const Vector& x);
	real log_pdf(const Vector& x, const Vector& y)
	{
		return - (y - expected_value(x)).SquareNorm();
	}
    void BootstrapBandwidth();
    void Show()
    {
    }
protected:
    int nearest_neighbour_size; ///< what size to use for the nearest neighbour
    KDTree<PointPair> kd_tree; ///< The tree, for faster access
    
    
};



#endif
