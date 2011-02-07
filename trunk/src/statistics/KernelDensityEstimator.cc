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
                                               real initial_bandwidth,
                                               int knn)
    : n(n_dimensions),
      b(initial_bandwidth),
      change_b(true),
      nearest_neighbour_size(knn),
      kd_tree(n_dimensions)
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
    if (nearest_neighbour_size > 0) {
        kd_tree.AddVectorObject(x, &points.back());
    }
}

real KernelDensityEstimator::log_pdf(const Vector& x)
{
    real C = - 0.5 * ((real) n) * log(2.0 * M_PI);
    
    // If no points are stored, use a standard normal density
    if (!points.size()) {
        real d = x.Norm(2.0);
        real r = C - 0.5 * d * d;
        printf ("! %f %f\n", r, exp(r));
        return r;
        
    }

    // otherwise, do the kernel estimate
    real ib2 = 1.0 / (b * b);
    real log_P = LOG_ZERO;
    if (nearest_neighbour_size == 0) {
        for (std::list<WeightedPoint>::iterator it = points.begin();
             it != points.end();
             ++it) {
            real d = SquareNorm(&x, &(it->x));
            real log_p_i = C - 0.5 * d * ib2;
            log_P = logAdd(log_P, log_p_i);
        }
    } else {
        int K = points.size();
        if (K > nearest_neighbour_size) {
            K = nearest_neighbour_size;
        }
        OrderedFixedList<KDNode> node_list = kd_tree.FindKNearestNeighbours(x, K);
        
        for (std::list<std::pair<real, KDNode*> >::iterator it = node_list.S.begin();
             it != node_list.S.end();
             ++it) {
            KDNode* node = it->second;
            WeightedPoint* p = kd_tree.getObject(node);
            real d = SquareNorm(&x, &(p->x));
            real log_p_i = C - 0.5 * d * ib2;
            log_P = logAdd(log_P, log_p_i);
        }
    }

    real N = (real) points.size(); 
    return log_P - log(b) - log(N);
}

/// Use bootstrapping to estimate the bandwidth
void KernelDensityEstimator::BootstrapBandwidth()
{
    KernelDensityEstimator kde(n, b, nearest_neighbour_size);
    std::list<WeightedPoint> test_data;
    fprintf(stderr, "Boostrapping bandiwdth\n");
    for (std::list<WeightedPoint>::iterator it = points.begin();
         it != points.end();
         ++it) {
        if (urandom() < 0.3) {
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
            fprintf (stderr, "b: %f -> %f (%f %f) # bandwidth reduced\n",
                     b, current_b,
                     log_p, current_log_p);
            b = current_b;
            log_p = current_log_p;
        } else {
            fprintf (stderr, "b: %f = %f (%f %f) # bandwidth found\n",
                     b, current_b,
                     log_p, current_log_p);
            break;
        }
        
        current_b *= 0.5;

    }
}


