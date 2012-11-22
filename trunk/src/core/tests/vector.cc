/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "Vector.h"
#include "EasyClock.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

int main(int argc, char** argv)
{
    
    int N = 100000;
    int iter=100;
    Vector x(N);
    Vector y(N);

    gsl_vector* gx = gsl_vector_alloc(N);
    gsl_vector* gy = gsl_vector_alloc(N);

    for (int i=0; i<x.Size(); ++i) {
        x(i) = i;
        y(i) = 2*i - 1;
        gsl_vector_set(gx, i, x(i));
        gsl_vector_set(gy, i, y(i));
    }

    {
        double total_time = 0.0;
        for (int k=0; k<10000; ++k) {
            double start_time = GetCPU();
            Vector z = x * y;
            double end_time = GetCPU();
            total_time += end_time - start_time;
        }
        logmsg("total time %f\n", total_time);
    }



    {
        double total_time = 0.0;
        gsl_vector* gz = gsl_vector_alloc(N);
        for (int k=0; k<10000; ++k) {
            double start_time = GetCPU();
            //gsl_vector_memcpy(gz, gx);
            //gsl_vector_mul(gz, gy);
            double end_time = GetCPU();
            total_time += end_time - start_time;
        }
        gsl_vector_free(gz);
        logmsg("total time %f\n", total_time);
    }

    gsl_vector_free(gx);
    gsl_vector_free(gy);
    return 0;
}

#endif
