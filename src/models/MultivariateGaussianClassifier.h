/* -*- Mode: c++;  -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef MULTIVARIATE_GAUSSIAN_CLASSIFIER_H
#define MULTIVARIATE_GAUSSIAN_CLASSIFIER_H

#include "NormalDistribution.h"
#include <vector>

/// A classifier using a multivariate gaussian for each class.
class MultivariateGaussianClassifier
{
public:
    const int n_inputs; ///< dimensionality  of observations
    const int n_classes; ///< number of classes
    std::vector<MultivariateNormalUnknownMeanPrecision*> class_distribution; ///< the distribution for each class
    Vector output;
    MultivariateGaussianClassifier(int n_inputs_, int n_classes_);
    ~MultivariateGaussianClassifier();
    int Classify(const Vector& x) 
    {
        return ArgMax(Output(x));
    }
    Vector& Output(const Vector& x);
    real Observe(const Vector& x, const int label);

};

#endif
