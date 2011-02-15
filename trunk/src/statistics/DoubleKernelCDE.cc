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

#include "DoubleKernelCDE.h"
#include "Random.h"


DoubleKernelCDE::DoubleKernelCDE(int n_x_dimensions,
                int n_y_dimensions,
                real initial_bandwidth)
:
    n_x(n_x_dimensions),
    n_y(n_y_dimensions),
    b_x(initial_bandwidth),
    b_y(initial_bandwidth)
{
}

real DoubleKernelCDE::Observe(const Vector& x, const Vector& y)
{
    real p = pdf(x, y);
    AddPoint(x, y);
    return p;
}

real DoubleKernelCDE::pdf(const Vector& x, const Vector& y)
{
    return exp(log_pdf(x, y));
}

/**  Compute the log pdf.

     For the i-th point, calculate
     K_y(y - y_i)

     
 */
real DoubleKernelCDE::log_pdf(const Vector& x, const Vector& y)
{
    real C = - 0.5 * ((real) n_y) * log(2.0 * M_PI);
    
    // If no points are stored, use a standard normal density
    if (!D.size()) {
        real d = y.Norm(2.0);
        real r = C - 0.5 * d * d;
        printf ("! %f %f\n", r, exp(r));
        return r;
    }
        // otherwise, do the kernel estimate
    real ib2 = 1.0 / (b_y * b_y);
    real log_Z = LOG_ZERO;
    real log_P = LOG_ZERO;
    for (std::vector<PointPair>::iterator it = D.begin();
         it != D.end();
         ++it) {
        real d_x = SquareNorm(&x, &(it->x));
        real d_y = SquareNorm(&y, &(it->y));
        real log_p_i = C - 0.5 * d_y * ib2;
        real log_p_c = -d_x;
        log_P = logAdd(log_P, log_p_c + log_p_i);
        log_Z = logAdd(log_Z, log_p_c);
    }
    return log_P - log_Z - log(b_y);
}

void DoubleKernelCDE::BootstrapBandwidth()
{
    DoubleKernelCDE kde(n_x, n_y, b_x);
    std::vector<PointPair> test_data;
    fprintf(stderr, "Bootstrapping bandiwdth\n");
    for (std::vector<PointPair>::iterator it = D.begin();
         it != D.end();
         ++it) {
        if (urandom() < 0.3) {
            test_data.push_back(*it);
        } else {
            kde.AddPoint(it->x, it->y);
        }
    }

    real current_b = b_x;
    real log_p = LOG_ZERO;
    while (1) {
        kde.b_x = current_b;
        kde.b_y = current_b;
        // Get log-likelihood of b.
        real current_log_p = 0;
        for (std::vector<PointPair>::iterator p = test_data.begin();
             p != test_data.end();
             ++p) {
            current_log_p += kde.log_pdf(p->x, p->y);
        }        

        if (current_log_p > log_p) {
            fprintf (stderr, "b: %f -> %f (%f %f) # bandwidth reduced\n",
                     b_x, current_b,
                     log_p, current_log_p);
            b_x = current_b;
            b_y = current_b;
            log_p = current_log_p;
        } else {
            fprintf (stderr, "b: %f = %f (%f %f) # bandwidth found\n",
                     b_x, current_b,
                     log_p, current_log_p);
            break;
        }
        
        current_b *= 0.5;

    }

    
}


