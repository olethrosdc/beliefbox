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
#include "KNNClassifier.h"
#include "MultivariateGaussianClassifier.h"
#include "ReadFile.h"
#include "Matrix.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>

template <class C>
void Evaluate(C& classifier, Matrix& data, std::vector<int>& labels)
{
    real n_errors = 0;
    real accuracy = 0;
    int T = data.Rows();
    for (int t=0; t<T; ++t) {
        Vector x = data.getRow(t);
        if (classifier.Classify(x) != labels[t]) {
            n_errors+=1;
        }
        accuracy += classifier.output(labels[t]);
    }
    printf ("%f %f\n", n_errors / (real) T, accuracy / (real) T);
}

int main(int argc, char** argv)
{
    Matrix data; 
    std::vector<int> labels;

    Matrix test_data; 
    std::vector<int> test_labels;
    
    printf("Arguments: training_data, n_classifiers, test_data\n");
    ReadClassData(data, labels, argv[1]);

    //fprintf(stderr, "Reading %s\n", argv[1]);
    int n_classifiers = atoi(argv[2]);
    bool test = false;
    if (argc==4) {
        //fprintf(stderr, "Reading %s\n", argv[3]);
        ReadClassData(test_data, test_labels, argv[3]);
        test = true;
    }
    int T = data.Rows();
    int n_inputs = data.Columns();
    int n_classes = 1 + Span(labels);

    Vector mean(n_inputs);
    Vector std(n_inputs);

    for (int i=0; i<n_inputs; ++i) {
        mean(i) = 0;
        std(i) = 0;
        for (int t=0; t<T; ++t) {
            mean(i) += data(t,i);
        }
        mean *= (1.0 / (real) T);
        for (int t=0; t<T; ++t) {
            real delta = data(t,i) - mean(i);
            std(i) += delta * delta;
        }
        std(i) *= (1.0 / (real) T);
    }
    
    
    for (int i=0; i<n_inputs; ++i) {
        real istd = 1;
        if (std(i) >= 0.001) {
            istd = 1 / std(i);
        }
        for (int t=0; t<data.Rows(); ++t) {
            data(t,i) -= mean(i);
            data(t,i) *= istd;                
        }
    }
    
    //printf ("Mean: "); mean.print(stdout);
    //printf ("Std:  "); std.print(stdout);
    //printf ("Matrix:\n");
    //data.print(stdout);


    if (test) {
        for (int i=0; i<n_inputs; ++i) {
            real istd = 1;
            if (std(i) >= 0.001) {
                istd = 1 / std(i);
            }
            for (int t=0; t<test_data.Rows(); ++t) {
                test_data(t,i) -= mean(i);
                test_data(t,i) *= istd;                
            }
        }
    }

    //LinearClassifier classifier(n_inputs, n_classes);
    //LinearClassifierMixture classifier(n_inputs, n_classes, n_classifiers);
    //HashedLinearClassifierMixture classifier(n_inputs, n_classes, n_classifiers);
    //MultivariateGaussianClassifier classifier(n_inputs, n_classes);
    KNNClassifier classifier(n_inputs, n_classes, 3);

    int n_iter = 1;
    real alpha = 0.001;

    //classifier.setStepSize(alpha);

    printf ("# K: %d, T: %d, d: %d inputs, n: %d classes, iter:%d, alpha: %f\n",
            n_classifiers,
            T,
            n_inputs,
            n_classes,
            n_iter,
            alpha);


    for (int iter=0; iter<n_iter; ++iter) {
        real n_errors = 0;
        real accuracy = 0;
        for (int t=0; t<T; ++t) {
            Vector x = data.getRow(t);
            ///Hash
            if (classifier.Classify(x) != labels[t]) {
                n_errors+=1;
            }
            accuracy += classifier.output(labels[t]);
            classifier.Observe(x, labels[t]);
            //printf ("%f %f\n", n_errors / (real) t, accuracy / (real) t);
        }
        //printf ("%f %f\n", n_errors / (real) T, accuracy / (real) T);
    }
    
    printf ("# TRAIN \n");
    if (1) {
        Evaluate(classifier, data, labels);
    }

    printf ("# TEST \n");
    if (test) {
        Evaluate(classifier, test_data, test_labels);
    }
}



#endif
