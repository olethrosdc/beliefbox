/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef MULTIVARIATE_NORMAL_H
#define MULTIVARIATE_NORMAL_H

#include "NormalDistribution.h"

/// Multivariate Gaussian probability distribution
class MultivariateNormal : public VectorDistribution
{
 private:
    int n_dim;
    Vector mean;
    Matrix accuracy;
    real determinant;
 public:
    MultivariateNormal(const int n_dim_);
    MultivariateNormal(const Vector& mean_, const Matrix& accuracy_);
    void setMean(const Vector& mean_)
    {
        mean = mean_;
    }
    void setAccuracy(const Matrix& accuracy_)
    {
        accuracy = accuracy_;
        accuracy.LUDecomposition(determinant);
    }
    virtual ~MultivariateNormal() {}
    virtual void generate(Vector& x) const;
    virtual Vector generate() const;
    virtual real log_pdf(const Vector& x) const;
    virtual real pdf(const Vector& x) const
    {
        return exp(log_pdf(x));
    }
    void Show() const;
};


class Student;


#endif
