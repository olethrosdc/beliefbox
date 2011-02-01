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

#include "KernelDensityEstimator.h"

/// Constructor
KernelDensityEstimator::KernelDensityEstimator(int n_dimensions,
                                               real initial_bandwidth)
    : n(n_dimensions),
      b(initial_bandwidth),
      change_b(true)
{
    
}


real KernelDensityEstimator::Observe(const Vector& x)
{
    real p = pdf(x);
    AddPoint(x);
    return p;
}


/// Add a point x with weight w (defaults to w = 1)
void KernelDensityEstimator::AddPoint(const Vector& x,  real w)
{
    points.push_back(WeightedPoint(x, w));
}

real KernelDensityEstimator::log_pdf(const Vector& x)
{
    real C = - 0.5 * ((real) n) * log(2.0 * M_PI);
    
    // If no points are stored, use a standard normal density
    if (!points.size()) {
        real d = x.Norm(2.0);
        return C - 0.5 * d * d;
        
    }

    // otherwise, do the kernel estimate
    real ib2 = 1.0 / (b * b);
    real log_P = LOG_ZERO;
    for (std::list<WeightedPoint>::iterator it = points.begin();
         it != points.end();
         ++it) {
        real d = SquareNorm(&x, &(it->x));
        real log_p_i = C - 0.5 * d * ib2;
        log_P = logAdd(log_P, log_p_i);
    }
    real N = (real) points.size();
    return log_P - log(b) - log(N);
}

/// Return the log pdf at x
real KernelDensityEstimator::pdf(const Vector& x)
{
    //printf ("%f %f\n", b, log_pdf(x));
    return exp(log_pdf(x));
}

/// Use bootstrapping to estimate the bandwidth
void KernelDensityEstimator::BootstrapBandwidth()
{
    KernelDensityEstimator kde(n, b);
    std::list<WeightedPoint> test_data;
    for (std::list<WeightedPoint>::iterator it = points.begin();
         it != points.end();
         ++it) {
        if (urandom() < 1.0/3.0) {
            test_data.push_back(*it);
        } else {
            kde.AddPoint(it->x);
        }
    }

    real current_b = b;
    real log_p = LOG_ZERO;
    while (1) {
        kde.b = current_b;

        // Get log-likelihood of b.
        real current_log_p = 0;
        for (std::list<WeightedPoint>::iterator p = test_data.begin();
             p != test_data.end();
             ++p) {
            current_log_p += kde.log_pdf(p->x);
            //current_log_p += kde.pdf(p->x);
        }        

        if (current_log_p > log_p) {
            fprintf (stderr, "b: %f -> %f (%f %f)\n",
                     b, current_b,
                     log_p, current_log_p);
            b = current_b;
            log_p = current_log_p;
        } else {
            fprintf (stderr, "b: %f = %f (%f %f)\n",
                     b, current_b,
                     log_p, current_log_p);
            break;
        }
        
        current_b *= 0.5;

    }
}
