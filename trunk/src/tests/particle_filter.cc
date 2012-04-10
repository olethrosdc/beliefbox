/* -*- Mode: C++; -*- */
/* VER: $Id: test_particle_filter.c,v 1.2 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifdef MAKE_MAIN

#include "ParticleFilter.h"
#include "ActionValueEstimate.h"
#include "SingularDistribution.h"
#include <cstdio>
#include <cstdlib>



#if 0
int main (void)
{
    int N = 128;

    NormalDistribution actual_transitions(0.0f, 0.1f);
    NormalDistribution actual_observations(0.0f, 1.0f);
    NormalDistribution transitions(0.0f, 0.1f);
    NormalDistribution observations(0.0f, 1.0f);
        
    ParticleFilter pf(N, &transitions, &observations);
    PopulationEstimate pe(1, N, 0.01);

    int T = 1000;
    FILE* f = fopen("pf_test","w");
    if (!f) {
        fprintf(stderr, "Could not open file\n");
        exit(-1);
    }
    real EX = drand48();
    for (int t=0; t<T; t++) {
        real X = EX + actual_observations.generate();
        pf.Observe(X);
        pe.Observe(0, X);
        fprintf(f, "%f %f %f %f %f\n", EX, pf.GetMean(), pe.GetMean(0), pf.GetVar(), pe.GetVar(0));
        EX = EX + actual_transitions.generate();
    }
    fclose(f);
    return 0;
}
#else

int main (void)
{
    int N = 16;
    real EX=urandom();
    SingularDistribution actual_transitions(0.0f);
    BernoulliDistribution actual_observations(EX);
    UniformDistribution transitions(-0.1f, 0.1f);
    BernoulliDistribution observations;
    UniformDistribution prior(0.0f, 1.0f);
                                            
        
    BernoulliGridParticleFilter pf(N, &prior, &transitions, &observations);

    int T = 1000;
    FILE* f = fopen("pf_test","w");
    FILE* values = fopen("pf_values","w");

    if (!f || !values) {
        fprintf(stderr, "Could not open file\n");
        exit(-1);
    }

    for (int t=0; t<T; t++) {
        real X = actual_observations.generate();
        pf.Observe(X);
        fprintf(f, "%f %f %f\n", EX, pf.GetMean(), pf.GetVar());
        for (int n=0; n<N; n++) {
            fprintf (values, "%f ", pf.y[n]);
        }
        for (int n=0; n<N; n++) {
            fprintf (values, "%f ", pf.w[n]);
        }
        fprintf (values, "\n");
    }
    fclose (values); 
    fclose(f);
    return 0;
}

#endif
#endif
