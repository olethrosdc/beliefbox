/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "LinearClassifier.h"
#include "ClassifierMixture.h"
#include "ReadFile.h"
#include "Matrix.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>

int main(int argc, char** argv)
{
    Matrix data; 
    std::vector<int> labels;
    
    ReadClassData(data, labels, argv[1]);
    
    int T = data.Rows();
    int n_inputs = data.Columns();
    int n_classes = 1 + Span(labels);
    //LinearClassifier classifier(n_inputs, n_classes);
    LinearClassifierMixture classifier(n_inputs, n_classes, 2);

    int n_iter = 128;
    real alpha = 0.001;
    classifier.setStepSize(alpha);
    for (int iter=0; iter<n_iter; ++iter) {
        real n_errors = 0;
        real accuracy = 0;
        for (int t=0; t<T; ++t) {
            Vector x = data.getRow(t);
            if (classifier.Classify(x) != labels[t]) {
                n_errors+=1;
            }
            accuracy += classifier.output(labels[t]);
            classifier.Observe(x, labels[t]);
        }
        printf ("%f %f)\n", n_errors / (real) T, accuracy / (real) T);
    }

    classifier.Show();

    
}

#endif
