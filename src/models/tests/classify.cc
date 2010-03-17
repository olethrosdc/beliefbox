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
    LinearClassifier classifier(n_inputs, n_classes);

    int n_iter = 1000;
    real alpha = 0.01;
    StochasticGradientClassifier sgc(classifier);
    sgc.setStepSize(alpha);
    for (int iter=0; iter<n_iter; ++iter) {
        real n_errors = 0;
        for (int t=0; t<T; ++t) {
            Vector x = data.getRow(t);
            if (classifier.Classify(x) != labels[t]) {
                n_errors+=1;
            }
            sgc.Observe(x, labels[t]);
        }
        printf ("%f (%d)\n", n_errors / (real) T, T);
    }
    classifier.Show();
}

#endif
