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

#include "LinearClassifier.h"
#include "Random.h"

LinearClassifier::LinearClassifier(int n_inputs_, int n_classes_)
    : n_inputs(n_inputs_), n_classes(n_classes_),
      params(n_inputs, n_classes), bias(n_classes),
      output(n_classes),
      dc_dg(n_classes), dg_df(n_classes),
      alpha(0.1)
{
    for (int j=0; j<n_classes; ++j) {
        for (int i=0; i<n_inputs; ++i) {
            params(i,j) = urandom() - 0.5;
        }
        bias(j) = urandom() - 0.5;
    }
}

Vector& LinearClassifier::Output(const Vector& x)
{
    Matrix tmp(n_classes, n_inputs);
    for (int i=0; i<n_inputs; ++i) {
        for (int j=0; j<n_classes; ++j) {
            tmp(j,i) = params(i,j);
        }
    }
    hidden = ((const Matrix&) tmp)* x + bias;
    SoftMax(hidden, output, 1.0);
    return output;
}

void LinearClassifier::Show()
{
    params.print(stdout);
}


void LinearClassifier::Observe(const Vector& x, int label)
{
    assert(x.Size() == n_inputs);
    assert(label >= 0 && label < n_classes);
    Output(x);
#if 0
    printf("# x: "); x.print(stdout);
    printf("# f: "); hidden.print(stdout);
    printf("# g: "); output.print(stdout);
#endif
    dc_dg.Clear();
    dc_dg(label) = 1;

    real sum = 0.0;
    for (int j=0; j<n_classes; ++j) {
        //dc_dg(j) -= output(j);
        sum += dc_dg(j);
    }

    for (int j=0; j<n_classes; ++j) {
        //        dg_df(j) = hidden(j) * output(j) * (1 + output(j));
        dg_df(j) = dc_dg(j) - output(j) * sum;
    }

    for (int j=0; j<n_classes; ++j) {
        real delta = alpha * dc_dg(j) * dg_df(j);
        for (int i=0; i<n_inputs; ++i) {
            params(i, j) += delta * x(i);
        }
        bias(j) += delta;
    }

}



