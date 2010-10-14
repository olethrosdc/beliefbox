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


#ifndef SPARSE_LINEAR_CLASSIFIER_H
#define SPARSE_LINEAR_CLASSIFIER_H

#include "real.h"
#include "Vector.h"
#include "Matrix.h"
#include "LinearClassifier.h"
#include "Random.h"

/** A sparse linear classifier.
    
    This is a sparse version of the linear classifier.
    It uses a \f$m \times n\f$ projection matrix \f$M\f$ to obtain
    \f[
    h = Mx,
    \f]
    with \f$x \in R^n\f$, \f$h \in R^m\f$.
 */
class SparseLinearClassifier
{
public:
    const int n_inputs;
    const int n_classes;
    const int projection_size;
    Matrix M; ///< projection matrix
    SparseLinearClassifier(int n_inputs_, int n_classes_, int projection_size_)
        : n_inputs(n_inputs_),
          n_classes(n_classes_),
          projection_size(projection_size_),
          M(projection_size, n_inputs),
          classifier(projection_size, n_classes),
          output(classifier.output)
    {
        for (int i=0; i<projection_size; ++i) {
            for (int j=0; j<n_inputs; ++j) {
                M(i,j) = urandom() - 0.5;
            }
        }
    }
    int Classify(const Vector& x)
    {
        return classifier.Classify(tanh((const Matrix&) M * (const Vector&) x));
    }
    Vector& Output(const Vector& x)
    {
        return classifier.Output(tanh((const Matrix&) M * (const Vector&) x));
    }
    void Observe(const Vector& x, int label)
    {
        classifier.Observe(tanh((const Matrix&) M * (const Vector&) x), label);
    }
    void Show()
    {
        classifier.Show();
    }
    void setStepSize(real step_size) {
        classifier.setStepSize(step_size);
    }
protected:
    LinearClassifier classifier;
public:
    Vector& output;
};



#endif
