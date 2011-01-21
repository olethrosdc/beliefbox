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

#include"MultivariateGaussianClassifier.h"

/// Construct the class
MultivariateGaussianClassifier::MultivariateGaussianClassifier(int n_inputs_, int n_classes_)
    : n_inputs(n_inputs_),
      n_classes(n_classes_),
      class_distribution(n_classes),
      prior(n_classes, 0.5),
      output(n_classes)
{
    Vector mu(n_inputs);
    real tau = 1.0;
    real alpha = 1.0;
    Matrix T(Matrix::Unity(n_inputs, n_inputs));
    
    for (int i=0; i<n_classes; ++i) {
        class_distribution[i] = new MultivariateNormalUnknownMeanPrecision(mu, tau, alpha, T);
    }

}

/// Destroy
MultivariateGaussianClassifier::~MultivariateGaussianClassifier()
{
    for (int i=0; i<n_classes; ++i) {
        delete class_distribution[i];
    }
}

/// Classify a vector x
Vector& MultivariateGaussianClassifier::Output(const Vector& x)
{
    for (int i=0; i<n_classes; ++i) {
        output(i) = class_distribution[i]->log_pdf(x);
        if (isnan(output(i))) {
            Serror("Output %d is nan", i);
            output.print(stdout);
            class_distribution[i]->Show();
        }
    }
    output += log(prior.GetMean());
    real S = output.logSum();
    output -= S;
    for (int i=0; i<n_classes; ++i) {
        output(i) = exp(output(i));
        if (isnan(output(i))) {
            Serror("Output %d is nan", i);
            output.print(stdout);
            class_distribution[i]->Show();
            exit(-1);
        }
    }
    return output;
}

real MultivariateGaussianClassifier::Observe(const Vector& x, const int label)
{
    assert(label >= 0 && label < n_classes);
    real probability = Output(x)(label);
    class_distribution[label]->Observe(x);
    prior.Observe(label);
    return probability;
}

Vector MultivariateGaussianClassifier::getClassMean(const int label) const
{
    return class_distribution[label]->getMean();
}
