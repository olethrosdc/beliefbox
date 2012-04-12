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

#ifdef MAKE_MAIN
#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "SingularDistribution.h"
#include "BetaDistribution.h"


void TestData(std::vector<real>& data,
                    ConjugatePrior& prior);

int main (int argc, char** argv)
{
    BetaDistribution beta;
    NormalUnknownMeanPrecision normal;
    UnknownSingularDistribution fixed;

    int N = 1024;

    std::vector<real> data(N);
    for (int i=0; i<N; ++i) {
        data[i] = 0.5;
    }
    printf ("Beta\n");
    TestData(data, beta);
    printf ("Normal\n");
    TestData(data, normal);
    printf ("Fixed\n");
    TestData(data, fixed);
};


void TestData(std::vector<real>& data,
              ConjugatePrior& prior)
{
    int N = data.size();
    real log_p = 0;
    for (int i=0; i<N; ++i) {
        log_p += log(prior.Observe(data[i]));
    }
    
    printf("%f # log likelihood\n", log_p);

    for (int i=0; i<8; ++i) {
        printf ("%f ", prior.generate());
    }
    printf ("# generated data\n");
}

#endif

