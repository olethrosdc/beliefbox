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
    Vector hidden;
    Vector output;
    Vector dc_dg;
    Vector dg_df;
    real alpha;
    LinearClassifier(int n_inputs_, int n_classes_);
    int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector& Output(const Vector& x);
    void Observe(const Vector& x, int label);
    void Show();
    void setStepSize(real step_size) {
        alpha = step_size;
    }
};



#endif
