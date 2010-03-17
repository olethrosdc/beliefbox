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


#ifndef LINEAR_CLASSIFIER_H
#define LINEAR_CLASSIFIER_H

#include "real.h"
#include "Vector.h"
#include "Matrix.h"

class LinearClassifier
{
public:
    const int n_inputs;
    const int n_classes;
    Matrix params;
    Vector bias;
    LinearClassifier(int n_inputs_, int n_classes_);
    int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector Output(const Vector& x);
    void Show();
};

class LinearClassifierMixture
{
public:
    const int n_inputs;
    const int n_classes;
    std::vector<LinearClassifier*> classifiers;
    Vector w;
    LinearClassifierMixture(std::vector<LinearClassifier*> classifiers_,
                            Vector w_)
    {
    }
    int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector Output(const Vector& x);
    void Show();
};

class StochasticGradientClassifier
{
protected:
    LinearClassifier& classifier;
    real alpha;
    Matrix& params;
    Vector& bias;
public:
    StochasticGradientClassifier(LinearClassifier& classifier_) : 
        classifier(classifier_),
        alpha(0.1),
        params(classifier.params),
        bias(classifier.bias)
    {
    }
    void setStepSize(real alpha_) {
        assert(alpha_ >= 0);
        alpha = alpha_;
    }
    void Observe(const Vector& x, int label);
    
};

#endif
