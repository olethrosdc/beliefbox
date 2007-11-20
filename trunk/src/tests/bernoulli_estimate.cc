/* -*- Mode: C++; -*- */
/* VER: $Id: test_bernoulli_estimate.c,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "DiscreteBanditPolicy.h"
#include "ActionValueEstimate.h"
#include "PFActionValueEstimate.h"
#include "real.h"
#include "Distribution.h"

void Evaluate(char* fname, Policy& policy, int n_actions);

// if the result is worse but it seems that it could be better then we have to
// calculate that.
int main (int argc, char** argv)
{

    if (argc!=4) {
        fprintf (stderr, "Args: n_iters n_samples alpha\n");
        exit(-1);
    }

    setRandomSeed(time(NULL));

    int n_iters = atoi(argv[1]);
    int n_samples = atoi(argv[2]);
    real alpha = atof(argv[3]);
    fprintf (stderr, "n_iters: %d n_samples: %d\n", 1, n_samples);
    PopulationEstimate population_estimate(1, n_samples, alpha);
    PopulationEstimate population_estimate2(1, n_samples, alpha*0.1f);
    PFActionValueEstimate pf_estimate(1, n_samples);
    BernoulliEstimate bernoulli_estimate (1, n_samples, 1.0f, 1.0f);

    population_estimate.Reset();
    population_estimate2.Reset();
    pf_estimate.Reset();
    bernoulli_estimate.Reset();
    
    real Ex = 0.57; //urandom();
    fprintf (stderr, "Ex: %f\n", Ex);
    printf ("# Ex: %f\n", Ex);
    
    for (int i=0; i<n_iters; i++) {
        real r = 0.0f;
        if (urandom() < Ex) {
            r = 1.0f;
        }
        population_estimate.Observe(0, r);
        population_estimate2.Observe(0, r);
        pf_estimate.Observe(0, r);
        bernoulli_estimate.Observe(0, r);
        printf ("%f %f %f %f %f\n",
                Ex,
                population_estimate.GetMean(0),
                population_estimate2.GetMean(0),
                pf_estimate.GetMean(0),
                bernoulli_estimate.GetMean(0));
    }
    return 0;
}
#endif
