/* -*- Mode: C++ -*- */
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

#include "BoundedLogNormalDistribution.h"
#include "BetaDistribution.h"
int main()
{
    int T = 1000;
    real a = 2;
    real b= 2;
    BetaDistribution beta(a, b);

    Vector U(1);
    Vector L(1);
    L(0) = 0;
    U(0) = 1;
    BoundedLogNormal bounded_log_normal(L, U);

    Vector x(1);
    
    // print out the transform
    for (real z = 0; z<=1; z+=0.01) {
        x(0) = z;
        Vector y = bounded_log_normal.transform(x);
        printf ("%f %f # z_y\n", z, y(0));
    }

    // print the prpobabilities
    for (int t=0; t<T; ++t) {
        x(0) = beta.generate();    
        real p = bounded_log_normal.Observe(x);
        printf ("%f # p_t\n", p);
    }

    // print out the transform
    for (real z = 0; z<=1; z+=0.01) {
        x(0) = z;
        real p = bounded_log_normal.pdf(x);
        printf ("%f %f # p_x\n", z, p);
    }


    
    
  
}
#endif
