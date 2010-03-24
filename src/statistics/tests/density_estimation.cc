/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "ContextTreeRealLine.h"
#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"


int main (int argc, char** argv)
{
    BetaDistribution distribution(100,200);
    ContextTreeRealLine pdf(2,0);
    int T = atoi(argv[1]);
    if (T < 0) {
        Serror("T should be >= 0\n");
        exit(-1);
    }
    for (int t=0; t<T; ++t) {
        //real x = urandom()*0.3 + 0.25;
        real x = distribution.generate();
        //std::cout << x << std::endl;
        real p = pdf.Observe(x);
        //        std::cout << p << std::endl;
    }
#if 1
    for (real x=0; x<1; x+=0.001) {
        real p = pdf.pdf(x);
        printf ("%f %f\n", x, p);
    }
    //pdf.Show();
#endif
    return 0;
}

#endif
