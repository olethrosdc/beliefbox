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
#include "SparseLinearClassifier.h"
#include "ClassifierMixture.h"
#include "KNNClassifier.h"
#include "MultivariateGaussianClassifier.h"
#include "ConditionalKDGaussianClassifier.h"
#include "ConditionalKDNNClassifier.h"
#include "MersenneTwister.h"
#include "ReadFile.h"
#include "Matrix.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <getopt.h>

template <class C>
void Evaluate(C& classifier, Matrix& data, std::vector<int>& labels, const char* s)
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
    printf ("%f %f #%s\n", n_errors / (real) T, accuracy / (real) T, s);
}

template <class C>
void Train(C& classifier, Matrix& data, std::vector<int>& labels, bool randomise=false)
{
    real n_errors = 0;
    real accuracy = 0;
    int T = data.Rows();
    std::vector<int> indices;
    if (randomise) {
        MersenneTwisterRNG rng;
        indices.resize(T);
        for (int i=0; i<T; ++i) {
            indices[i] = i;
        }
        for (int i=0; i<T; ++i) {
            int j = i + rng.discrete_uniform(T - i);
            int tmp = indices[j];
            indices[j] = indices[i];
            indices[i] = tmp;
        }
    }

    for (int i=0; i<T; ++i) {
        int t;
        if (randomise) {
            t = indices[i];
        } else {
            t = i;
        }
        Vector x = data.getRow(t);
        ///Hash
        if (classifier.Classify(x) != labels[t]) {
            n_errors+=1;
        }
        accuracy += classifier.output(labels[t]);
        classifier.Observe(x, labels[t]);
        //printf ("%f %f #PREQ\n", n_errors / (real) (T), accuracy / (real) (T));
    }
    printf ("%f %f #PREQ\n", n_errors / (real) (T), accuracy / (real) (T));
	//classifier.Show();
}

enum ClassifierType {
    GAUSSIAN_TREE,
    KNN_TREE,
    NEAREST_NEIGHBOUR,
    GAUSSIAN,
    HASHED_MIXTURE,
    LINEAR,
    SPARSE_LINEAR
};
static const char* const help_text = "Usage: ... \n\
 --train file        File to use for training\n\
 --test file         File to use for evaluation\n\
 --sparse N          Sparse linear classifier with projection size N\n\
 --tree N            Tree Gaussian classifier\n\
 --mixture N         Hashed mixture of N classifiers (default 1-NN)\n\
 --knn K             Use a K-nearest neighbour classifier as base\n\
 --gaussian          Use a Gaussian classifier as base\n\
 --iterations N      Number of iterations N on the training set\n\
 --step_size u       Step size u in [0,1] to use for some methods\n\
 --normalise         Make the data zero mean, unit variance\n\
 --randomise         Shuffle the data before presentation\n";

int main(int argc, char** argv)
{
    Matrix data; 
    std::vector<int> labels;

    Matrix test_data; 
    std::vector<int> test_labels;
    
    printf("Arguments: training_data, n_classifiers, test_data\n");

    const char* train_filename = NULL;
    const char* test_filename = NULL;
    int tree_depth = 1;
    int n_classifiers = 1;
    int n_neighbours = 1;
    int n_iterations = 1;
    real step_size = 0.01;
    bool normalise = false;
    bool randomise = false;
    int projection_size = 4;
    bool use_knn = false;
    bool use_gaussian = false;
    ClassifierType classifier_type = LINEAR;
    {
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"train", required_argument, 0, 0}, //0
                {"test", required_argument, 0, 0}, //1
                {"tree", required_argument, 0, 0}, //2
                {"mixture", required_argument, 0, 0}, //3
                {"knn", required_argument, 0, 0}, //4
                {"iterations", required_argument, 0, 0}, //5
                {"step_size", required_argument, 0, 0}, // 6
                {"normalise", no_argument, 0, 0}, // 7 
                {"randomise", no_argument, 0, 0}, // 8
                {"sparse", required_argument, 0, 0}, // 9
                {"gaussian", no_argument, 0, 0}, // 10
                {0, 0, 0, 0}
            };
            c = getopt_long (argc, argv, "",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
            case 0:
#if 0
                printf ("option %s (%d)", long_options[option_index].name, option_index);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
#endif
                switch (option_index) {
                case 0: train_filename = optarg; break;
                case 1: test_filename = optarg; break;
                case 2:
                    tree_depth = atoi(optarg); 
                    classifier_type = GAUSSIAN_TREE;
                    break;
                case 3:
                    n_classifiers = atoi(optarg); 
                    classifier_type = HASHED_MIXTURE;
                    break;
                case 4: 
                    n_neighbours = atoi(optarg); 
                    use_knn = true;
                    break;
                case 5:
                    n_iterations = atoi(optarg); 
                    break;
                case 6:
                    step_size = atof(optarg); 
                    break;
                case 7:
                    normalise = true;
                    break;
                case 8:
                    randomise = true;
                    break;
                case 9:
                    projection_size = atoi(optarg);
                    classifier_type = SPARSE_LINEAR;
                    break;
                case 10:
                    use_gaussian = true;
                    break;
                default:
                    fprintf (stderr, "%s", help_text);
                    exit(0);
                    break;
                }
                break;
            case '0':
            case '1':
            case '2':
                if (digit_optind != 0 && digit_optind != this_option_optind)
                    printf ("digits occur in two different argv-elements.\n");
            digit_optind = this_option_optind;
            printf ("option %c\n", c);
            break;
            default:
                std::cout << help_text;
                exit (-1);
            }
        }
        
        if (optind < argc) {
            printf ("non-option ARGV-elements: ");
            while (optind < argc) {
                printf ("%s ", argv[optind++]);
                
            }
            printf ("\n");
        }
    }

    
    if (use_gaussian && use_knn)  {
        Serror("Must specify at most one of {use_gaussian, use_knn}\n");
        exit(-1);
    }

    if (use_knn && classifier_type == GAUSSIAN_TREE) {
        classifier_type = KNN_TREE;
        logmsg("Setting classifier to KNN tree");
    }

    if (classifier_type == LINEAR) {
        if (use_gaussian) {
            classifier_type = GAUSSIAN;
            logmsg("Setting classifier to multivariate Gaussian");
        }
        if (use_knn) {
            classifier_type = NEAREST_NEIGHBOUR;
            logmsg("Setting classifier to k-NN");
        }
    }

    if (!ReadClassData(data, labels, train_filename)) {
        Serror("Could not read train data\n");
    }

    
    
    bool test = false;
    if (test_filename) {
        if (ReadClassData(test_data, test_labels, test_filename)) {
            test = true;
        } else {
            Serror("Could not read test data\n");
        }
    }
    int T = data.Rows();
    int n_inputs = data.Columns();
    int n_classes = 1 + Span(labels);


    Vector mean(n_inputs);
    Vector std(n_inputs);
    Vector lower_bound(n_inputs);
    Vector upper_bound(n_inputs);


    for (int i=0; i<n_inputs; ++i) {
        mean(i) = 0;
        std(i) = 0;
        lower_bound(i) = -1;
        upper_bound(i) = +1;
        for (int t=0; t<T; ++t) {
            mean(i) += data(t,i);
            lower_bound(i) = std::min(lower_bound(i), data(t, i));
            upper_bound(i) = std::max(upper_bound(i), data(t, i));
        }
        mean *= (1.0 / (real) T);
        if (normalise) {
            for (int t=0; t<T; ++t) {
                real delta = data(t,i) - mean(i);
                std(i) += delta * delta;
            }
        }
        std(i) *= (1.0 / (real) T);
    }
    
    
    if (normalise) {
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
    }
    
    //printf ("Mean: "); mean.print(stdout);
    //printf ("Std:  "); std.print(stdout);
    //printf ("Matrix:\n");
    //data.print(stdout);


    if (test && normalise) {
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


    LinearClassifier linear_classifier(n_inputs, n_classes);
    linear_classifier.setStepSize(step_size);
    
    SparseLinearClassifier sparse_linear_classifier(n_inputs, n_classes, projection_size);
    sparse_linear_classifier.setStepSize(step_size);
    //LinearClassifierMixture classifier(n_inputs, n_classes, n_classifiers);
    //HashedLinearClassifierMixture linear_classifier_mixture(n_inputs, n_classes, n_classifiers, projection_size);
	std::vector<KNNClassifier*> experts(n_classifiers);
	for (int i=0; i<n_classifiers; ++i) {
		experts[i] = new KNNClassifier(n_inputs, n_classes, n_neighbours);
	}
    HashedClassifierMixture<KNNClassifier> hashed_classifier_mixture(n_inputs, n_classes, experts);

    MultivariateGaussianClassifier gaussian_classifier(n_inputs, n_classes);
    KNNClassifier knn_classifier(n_inputs, n_classes, n_neighbours);
    ConditionalKDGaussianClassifier tree_gaussian(2, tree_depth, lower_bound, upper_bound, n_classes);
    ConditionalKDNNClassifier tree_knn(2, tree_depth, lower_bound, upper_bound, n_classes);


    //classifier.setStepSize(alpha);

    printf ("# K: %d, T: %d, d: %d inputs, n: %d classes\n",
            n_classifiers,
            T,
            n_inputs,
            n_classes);
    
    for (int iter=0; iter<n_iterations; ++iter) {
        switch(classifier_type) {
        case GAUSSIAN_TREE:
            Train(tree_gaussian, data, labels, randomise);
            break;
        case KNN_TREE:
            Train(tree_knn, data, labels, randomise);
            break;
        case GAUSSIAN:
            Train(gaussian_classifier, data, labels, randomise);
            break;
        case NEAREST_NEIGHBOUR:
            Train(knn_classifier, data, labels, randomise);
            break;
        case HASHED_MIXTURE:
            Train(hashed_classifier_mixture, data, labels, randomise);
            break;
        case LINEAR:
            Train(linear_classifier, data, labels, randomise);
            break;
        case SPARSE_LINEAR:
            Train(sparse_linear_classifier, data, labels, randomise);
            break;
        }
    }
	

    
    if (1) {
        printf("# Evaluating ...\n");
        switch(classifier_type) {
        case KNN_TREE:
            Evaluate(tree_knn, data, labels, "TRAIN");
            break;
        case GAUSSIAN:
            Evaluate(gaussian_classifier, data, labels, "TRAIN");
            break;
        case GAUSSIAN_TREE:
            Evaluate(tree_gaussian, data, labels, "TRAIN");
            break;
        case NEAREST_NEIGHBOUR:
            Evaluate(knn_classifier, data, labels, "TRAIN");
            break;
        case HASHED_MIXTURE:
            Evaluate(hashed_classifier_mixture, data, labels, "TRAIN");
            break;
        case LINEAR:
            Evaluate(linear_classifier, data, labels, "TRAIN");
            break;
        case SPARSE_LINEAR:
            Evaluate(sparse_linear_classifier, data, labels, "TRAIN");
            break;
        }
    }

    if (test) {
        printf("# Evaluating ...\n");
        switch(classifier_type) {
        case KNN_TREE:
            Evaluate(tree_knn, test_data, test_labels, "TEST");
            break;
        case GAUSSIAN:
            Evaluate(gaussian_classifier, test_data, test_labels, "TEST");
            break;
        case GAUSSIAN_TREE:
            Evaluate(tree_gaussian, test_data, test_labels, "TEST");
            break;
        case NEAREST_NEIGHBOUR:
            Evaluate(knn_classifier, test_data, test_labels, "TEST");
            break;
        case HASHED_MIXTURE:
            Evaluate(hashed_classifier_mixture, test_data, test_labels, "TEST");
            break;
        case LINEAR:
            Evaluate(linear_classifier, test_data, test_labels, "TEST");
            break;
        case SPARSE_LINEAR:
            Evaluate(sparse_linear_classifier, test_data, test_labels, "TEST");
            break;
        }

    }

	for (int i=0; i<n_classifiers; ++i) {
		delete experts[i];
	}
	return 0;
}



#endif
