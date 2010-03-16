// -*- Mode: c++ -*-

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
    int Classify(const Vector x)
    {
        return ArgMax(Output(x));
    }
    Vector Output(const Vector x);
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
    void Observe(Vector x, int label);

};

#endif
