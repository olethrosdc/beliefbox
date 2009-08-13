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

#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <vector>
#include <cassert>
#include "Vector.h"

class RBF
{
public:
    Vector center;
    real beta;
    RBF(Vector& c, real b)
    {
        center.Resize(c.Size());
        center = c;
        beta = b;
        assert(b > 0);
    }
    real Evaluate(Vector& x)
    {
        real d = EuclideanNorm(&x, &center);
        return exp(-beta*d);
    }
    real logEvaluate(Vector& x)
    {
        return - beta * EuclideanNorm(&x, &center);
    }
};

class RBFBasisSet
{
protected:
    std::vector<RBF*> centers;
    std::vector<real> log_features;
    std::vector<real> features;
    bool valid_features;
    bool valid_log_features;
    int n_bases;
public:
    RBFBasisSet()
    {
        valid_features = false;
        valid_log_features = false;
        n_bases = 0;
    }
    void AddCenter(Vector& v, real b);
    void Evaluate(Vector& x);
    void logEvaluate(Vector& x);
    int size()
    {
        return n_bases;
    }
    real log_F(int j)
    {
        assert(j>=0 && j < n_bases);
        assert(valid_log_features);
        return log_features[j];
    }
    real F(int j)
    {
        assert(j>=0 && j < n_bases);
        assert(valid_features);
        return features[j];
    }

};
#endif
