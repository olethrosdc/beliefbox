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

#include "BasisSet.h"

RBFBasisSet::RBFBasisSet(const EvenGrid& grid, real bandwidth)
{
    for (int i=0; i<grid.getNIntervals(); ++i) {
        AddCenter(grid.getCenter(i), bandwidth);
    }
}

void RBFBasisSet::AddCenter(const Vector& v, real b)
{
    RBF* rbf = new RBF(v, b);
    centers.push_back(rbf);
    features.push_back(0.0);
    log_features.push_back(0.0);
    valid_features = false;
    valid_log_features = false;
    n_bases++;
}

void RBFBasisSet::logEvaluate(const Vector& x)
{
    real log_sum = LOG_ZERO;
    for (int i=0; i<n_bases; ++i) {
        log_features[i] = centers[i]->logEvaluate(x);
        log_sum = logAdd(log_features[i], log_sum);
    }
    for (int i=0; i<n_bases; ++i) {
        log_features[i] -= log_sum;
    }
    valid_log_features = true;
    valid_features = false;
}

void RBFBasisSet::Evaluate(const Vector& x)
{
    logEvaluate(x);
    for (int i=0; i<n_bases; ++i) {
        features[i] = exp(log_features[i]);
    }
    valid_features = true;
}

