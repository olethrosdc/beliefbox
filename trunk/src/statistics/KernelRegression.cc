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

#include "KernelRegression.h"

/// Constructor
KernelRegression::KernelRegression(int n_dimensions_x,
								   int n_dimensions_y,
								   real initial_bandwidth,
								   int knn)
    : n_x(n_dimensions_x),
	  n_y(n_dimensions_y),
      b(initial_bandwidth),
      change_b(true),
      nearest_neighbour_size(knn),
      kd_tree(n_dimensions_x)
{
    
}


real KernelRegression::Observe(const Vector& x, const Vector& y)
{
    //real p = pdf(x);
    AddPoint(x, y);
    return 1.0;
}


/// Add a point x with weight w (defaults to w = 1)
void KernelRegression::AddPoint(const Vector& x,  const Vector& y)
{
    points.push_back(PointPair(x, y));
    if (nearest_neighbour_size > 0) {
        kd_tree.AddVectorObject(x, &points.back());
    }
}

Vector KernelRegression::expected_value(const Vector& x) 
{

    Vector y(n_y);

    // If no points are stored, use the zero vector
    if (!points.size()) {
        return y;
    }

    // otherwise, do the kernel estimate
	real log_P = LOG_ZERO;

    real ib2 = 1.0 / (b * b);
    if (nearest_neighbour_size == 0) {
        for (std::list<PointPair>::const_iterator it = points.begin();
             it != points.end();
             ++it) {
            real d = SquareNorm(&x, &(it->x));
            real log_p_i = - 0.5 * d * ib2;
            log_P = logAdd(log_P, log_p_i);
			real p_i = exp(log_p_i);
			y += it->y * p_i;
		}
    } else {
        int K = points.size();
        if (K > nearest_neighbour_size) {
            K = nearest_neighbour_size;
        }
        OrderedFixedList<KDNode> node_list = kd_tree.FindKNearestNeighbours(x, K);
		real max_log_p_i;
		{
			std::list<std::pair<real, KDNode*> >::iterator it = node_list.S.begin();
			KDNode* node = it->second;
            PointPair* p = kd_tree.getObject(node);
			real d = SquareNorm(&x, &(p->x));
			max_log_p_i =  - 0.5 * d * ib2;
		}
        for (std::list<std::pair<real, KDNode*> >::iterator it = node_list.S.begin();
             it != node_list.S.end();
             ++it) {
            KDNode* node = it->second;
            PointPair* p = kd_tree.getObject(node);
            real d = SquareNorm(&x, &(p->x));
            real log_p_i =  - 0.5 * d * ib2 - max_log_p_i;
			real p_i = exp(log_p_i);
            log_P = logAdd(log_P, log_p_i);
			y += p->y * p_i;
        }
    }
	return y * exp(-log_P);
	//    return exp(log_Y - log_P);
}

/// Use bootstrapping to estimate the bandwidth
void KernelRegression::BootstrapBandwidth()
{
    KernelRegression kde(n_x, n_y, b, nearest_neighbour_size);
    std::list<PointPair> test_data;
    fprintf(stderr, "Boostrapping bandiwdth\n");
    for (std::list<PointPair>::iterator it = points.begin();
         it != points.end();
         ++it) {
        if (urandom() < 0.3) {
            test_data.push_back(*it);
        } else {
            kde.AddPoint(it->x, it->y);
        }
    }

    real current_b = b;
    real log_p = LOG_ZERO;
    while (1) {
        kde.b = current_b;
        // Get log-likelihood of b.
        real current_log_p = 0;
        for (std::list<PointPair>::iterator p = test_data.begin();
             p != test_data.end();
             ++p) {
            current_log_p += kde.log_pdf(p->x, p->y);
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


