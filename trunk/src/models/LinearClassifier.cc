#include "LinearClassifier.h"

LinearClassifier::LinearClassifier(int n_inputs_, int n_classes_)
    : n_inputs(n_inputs_), n_classes(n_classes_),
      params(n_inputs, n_classes), bias(n_classes)
{
    
}

Vector LinearClassifier::Output(const Vector x)
{
    const Matrix& w = params;
    return w * x + bias;
}


void StochasticGradientClassifier::Observe(Vector x, int label)
{
    assert(x.Size() == classifier.n_inputs);
    assert(label >= 0 && label < classifier.n_classes);

    for (int i=0; i<classifier.n_inputs; ++i) {
        params(i, label) += alpha * x(i);
    }
    bias(label) += alpha;
}

