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

#include "ClassifierMixture.h"
#include "Distribution.h"
#include "Random.h"

LinearClassifierMixture::LinearClassifierMixture(int n_inputs_,
                                                 int n_classes_,
                                                 int n_classifiers) 
    : n_inputs(n_inputs_),
      n_classes(n_classes_),
      classifiers(n_classifiers),
      w(n_classifiers),
      P(n_classifiers),
      output(n_classes),
      alpha(0.1)
{
    real p = 1.0 / (real) n_classifiers;
    for (int i=0; i<n_classifiers; ++i) {
        classifiers[i] = new LinearClassifier(n_inputs, n_classes);
        classifiers[i]->setStepSize(alpha);
        w[i] = p;
    }

}
LinearClassifierMixture::~LinearClassifierMixture()
{
    for (int i=0; i<(int) classifiers.size(); ++i) {
        delete classifiers[i];
    }
}

Vector& LinearClassifierMixture::Output(const Vector& x)
{
    output.Clear();
    for (int i=0; i<(int) classifiers.size(); ++i) {
        output += classifiers[i]->Output(x) * w[i];
    }
    return output;
}


void LinearClassifierMixture::Observe(const Vector& x, int label)
{
    assert(x.Size() == n_inputs);
    assert(label >= 0 && label < n_classes);
    real sum = 0;
    for (int i=0; i<(int) classifiers.size(); ++i) {
        P(i) = 0.1 + w(i) * (classifiers[i]->Output(x))(label);
        sum += P(i);
    }
    sum = 1.0 / sum;
    P *= sum;
    classifiers[DiscreteDistribution::generate(P)]->Observe(x, label);
}

void LinearClassifierMixture::setStepSize(real step_size)
{
    alpha = step_size;
    for (int i=0; i<(int) classifiers.size(); ++i) {
        classifiers[i]->setStepSize(alpha);
    }
    //printf ("# Setting step size to %f\n", alpha);
}

void LinearClassifierMixture::Show()
{
    for (int i=0; i<(int) classifiers.size(); ++i) {
        printf ("# %d\n", i);
        classifiers[i]->Show();
    }
    printf ("w: ");
    w.print(stdout);
}


HashedClassifierMixture::HashedClassifierMixture(int n_inputs_, int n_classes_, int n_classifiers) : LinearClassifierMixture(n_inputs_, n_classes_, n_classifiers)
{
    secret = true_random_bits(false);
}


void HashedClassifierMixture::Observe(const Vector& x, int label)
{
    assert(x.Size() == n_inputs);
    assert(label >= 0 && label < n_classes);

    unsigned long hash = secret;
    for (int i=0; i<x.Size(); ++i) {
        void* ptr = (void*) &x[i];
        unsigned long xi = *((unsigned long*) ptr);
        hash ^= xi;
    }

    for (int i=0; i<(int) classifiers.size(); ++i) {
        if (hash & 11) {
            //printf ("1");
            classifiers[i]->Observe(x, label);
        } else {
            //printf("0");
        }
        hash = hash >> 1;
    }
}
