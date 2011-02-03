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

#include "KernelConditionalDensityEstimator.h"

KernelConditionalDensityEstimator::KernelConditionalDensityEstimator(int n_x_dimensions,
                                                                     int n_y_dimensions,
                                                                     
                                                                     real initial_bandwidth)
    : n_x(n_x_dimensions),
      n_y(n_y_dimensions),
      p_xy(n_x + n_y, initial_bandwidth),
      p_x(n_x, initial_bandwidth)
{
}

real KernelConditionalDensityEstimator::Observe(const Vector& x, const Vector& y)
{
    real P_xy = p_xy.Observe(Concatenate(x, y));
    real P_x = p_x.Observe(x);
    return P_xy / P_x;
}

real KernelConditionalDensityEstimator::pdf(const Vector& x, const Vector& y)
{
    return exp(log_pdf(x, y));
}


real KernelConditionalDensityEstimator::log_pdf(const Vector& x, const Vector& y)
{
    real log_P_xy = p_xy.log_pdf(Concatenate(x, y));
    real log_P_x = p_x.log_pdf(x);
    return log_P_xy - log_P_x;
}
