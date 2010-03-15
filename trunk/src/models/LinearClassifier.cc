#include "LinearClassifier.h"

LinearClassifier::LinearClassifier(int n_inputs_, int n_classes_)
    : n_inputs(n_inputs_), n_classes(n_classes_), params(n_inputs, n_classes)
{
    
}

int LinearClassifier::classify(const Vector x)
{
    const Matrix& w = params;
    Vector y = w * x;
    return ArgMax(y);
}
