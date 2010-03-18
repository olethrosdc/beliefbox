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

#ifndef CLASSIFIER_MIXTURE_H
#define CLASSIFIER_MIXTURE_H

#include "LinearClassifier.h"
#include <vector>

class LinearClassifierMixture
{
public:
    const int n_inputs;
    const int n_classes;
    std::vector<LinearClassifier*> classifiers;
    Vector w; ///< classifier weights
    Vector P; ///< classifier selection probabilities
    Vector output;
    real alpha;
    LinearClassifierMixture(int n_inputs_, int n_classes_,
                            int n_classifiers);
    ~LinearClassifierMixture();
    int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector& Output(const Vector& x);
    void Observe(const Vector& x, int label);
    void Show();
    void setStepSize(real step_size);
};

#endif
