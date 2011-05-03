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

#ifndef MEAN_ESTIMATOR_H
#define MEAN_ESTIMATOR_H

#include "real.h"

class MeanEstimator
{
 public:
    real mean;
    int n_samples;
    real min;
    real max;
    MeanEstimator()
        : mean(0.0), n_samples(0), min(0), max(0)
    {
    }
    MeanEstimator(real mean_, int n_samples_)
        : mean(mean_), n_samples(n_samples_), min(mean_), max(mean_)
    {
    }
    void Reset(real mean_, int n_samples_)
    {
        mean = mean_;
        n_samples = n_samples_;
        min = mean;
        max = mean;
    }
    real Observe(real x)
    {
        real alpha = 1.0 / ((real) (++n_samples));
        mean = alpha * x + (1 - alpha)*mean;
        if (x < min) {
            min = x;
        } else if (x > max) {
            max = x;
        }
        return mean;
    }
    real GetMean() const
    {
        return mean;
    }
    real GenerateMarginal() const
    {
		real a = min;
		real b = max;
		if (urandom() < 0.5) {
			b = mean;
		} else {
			a = mean;
		}
		real u = urandom();
        return u * a + (1 - u) * b;
    }
    real Generate() const
    {
		real delta = urandom();
		real cb = (max - min) * sqrt(log(delta) / (2.0 * (real) n_samples));
		real u = urandom();
		return (mean - cb) * u + (mean + cb) * (1 - u);
    }
};





#endif
