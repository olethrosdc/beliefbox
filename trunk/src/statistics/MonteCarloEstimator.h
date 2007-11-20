/* -*- Mode: C++; -*- */
/* VER: $Id: MonteCarloEstimator.h,v 1.2 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MONTECARLO_ESTIMATOR_H
#define MONTECARLO_ESTIMATOR_H

#include "Sampling.h"
#include "Distribution.h"
#include <vector>

class MonteCarloEstimator
{
 public:
    int N;
    std::vector<real> y;
    std::vector<real> y2;
    std::vector<real> w;
    std::vector<real> w2;
    Distribution* transitions; ///< Transitions
    Distribution* observations; ///< Observations
	
    /// Constructor
    MonteCarloEstimator(int N, Distribution* prior, Distribution* T, Distribution* O)
    {
        this->N = N;
        y.resize(N);
        y2.resize(N);
        w.resize(N);
        w2.resize(N);
        real a = 1.0f/(real) N;
        for (int i=0; i<N; i++) {
            w[i] = a;
            y[i] = prior->generate();
        }
        this->transitions = T;
        this->observations = O;
    }
    
    void Observe(real x)
    {
        // Generate a set of samples from our current belief
        for (int n=0; n<N; n++) {
            int Yn = PropSample (w);
            y2[n] = y[Yn] + transitions->generate(); 
        }

        // Evaluate the new posterior up to a normalising constant
        real sum = 0.0f;
        for (int i=0; i<N; i++) {
            real likelihood = observations->pdf(x - y2[i]);
            real prior = 0.0f;
            for (int j=0; j<N; j++) {
                prior += w[j]*transitions->pdf(y[j] - y2[i]);
            }
            w2[i] = likelihood*prior;
            sum += w2[i];
        }

        // Normalise to create a new filtering distribution
        real isum = 1.0f / sum;
        for (int i=0; i<N; i++) {
            w[i] = w2[i] * isum;
            y[i] = y2[i];
        }
    }

    /// Get the current mean;
    real GetMean()
    {
        real mean = 0.0f;
        for (int i=0; i<N; i++) {
            mean += w[i] * y[i];
        }
        return mean;
    }
    
    /// Get the current mean;
    real GetVar()
    {
        real mean = GetMean();
        real var = 0.0f;
        for (int i=0; i<N; i++) {
            real delta = y[i] - mean;
            var += w[i] * delta * delta;
        }
        return var;
    }
};

#endif
