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

RBFBasisSet::RBFBasisSet(const EvenGrid& grid, real scale)
{
    n_bases = 0;
    for (int i=0; i<grid.getNIntervals(); ++i) {
        AddCenter(grid.getCenter(i), grid.delta * scale);
    }
    //logmsg("Added %d RBFs\n", centers.size());
}

RBFBasisSet::~RBFBasisSet()
{
    for (uint i=0;  i<centers.size(); ++i) {
        delete centers[i];
    }
}

void RBFBasisSet::AddCenter(const Vector& v, const Vector& b)
{
    RBF* rbf = new RBF(v, b);
    centers.push_back(rbf);
    n_bases++;
    features.Resize(centers.size());
    features[n_bases-1] = 0.0;
    log_features.Resize(centers.size());
    log_features[n_bases-1] = 0.0;
    valid_features = false;
    valid_log_features = false;
}

void RBFBasisSet::AddCenter(const Vector& v, real b)
{
    RBF* rbf = new RBF(v, b);
    centers.push_back(rbf);
    n_bases++;
    features.Resize(centers.size());
    features[n_bases-1] = 0.0;
    log_features.Resize(centers.size());
    log_features[n_bases-1] = 0.0;
    valid_features = false;
    valid_log_features = false;
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

//void RBFBasisSet::Evaluate(const Vector& x)
//{
//    logEvaluate(x);
//    for (int i=0; i<n_bases; ++i) {
//        features[i] = exp(log_features[i]);
//    }
//    valid_features = true;
//}

void RBFBasisSet::Evaluate(const Vector& x)
{
    for(int i = 0; i<n_bases; ++i){
        features[i] = centers[i]->Evaluate(x);
    }
    valid_log_features = false;
    valid_features = true;
}

