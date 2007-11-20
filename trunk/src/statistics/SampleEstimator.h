/* -*- Mode: C++; -*- */
/* VER: $Id: SampleEstimator.h,v 1.2 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SAMPLE_ESTIMATOR_H
#define SAMPLE_ESTIMATOR_H

#include <vector>
#include <cmath>
#include <cstdlib>

#include "real.h"

/** Sample estimate of a random variable.
    
A sample estimate of the mean of a random variable.
We use N samples to represent our current knowledge of the mean,
distributed according to some prior distribution.
 */
class SampleEstimator
{
public:
    std::vector<real> x; // expected reward of actions
    int N;
    real alpha;
    // Make an estimator
    SampleEstimator() 
    {
        alpha = 0.1f;
        SetNumberOfEstimates(16);
    }
    
    // Set the learning rate
    void SetLearningRate(real alpha)
    {
        this->alpha = alpha;
    }

    // randomly set estimates according to our prior
    void SetNumberOfEstimates(int N)
    {
        this->N = N;
        x.resize(N);
        Reset();
    }

    /// Reset 
    void Reset()
    {
        for (int i=0; i<N; i++) {
            x[i] = drand48();
        }
    }

    inline real Sample()
    {
        return x[rand()%N];
    }

    real GetMean()
    {
        real mean = 0.0f;
        for (int i=0; i<N; i++) {
            mean += x[i];
        }
        return mean / (real) N;
    }


    real GetVar()
    {
		real mean = GetMean();
		real var = 0.0f;
		for (int i=0; i<N; i++) {
			real delta = x[i] - mean;
			var += delta * delta;
		}
		return var / (real) N;
    }

    inline void Observe(real X)
    {
        for (int i=0; i<N; i++) {
            x[i] += alpha * (X - x[i]);
        }
    }
};

#endif
