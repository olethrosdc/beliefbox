// -*- Mode: c++ -*-

#ifndef LINEAR_CLASSIFIER_H
#define LINEAR_CLASSIFIER_H

#include "real.h"
#include "Vector.h"
#include "Matrix.h"

class LinearClassifier
{
protected:
    const int n_inputs;
    const int n_classes;
    Matrix params;
public:
    LinearClassifier(int n_inputs_, int n_classes_);
    int classify(const Vector x);
};

class StochasticTransient
{
protected:
    LinearClassifier& classifier;
public:
    
};

#endif
