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

#include "SparseLinearClassifier.h"
#include "Random.h"
#include <vector>

class LinearClassifierMixture 
{
public:
    const int n_inputs;
    const int n_classes;
    std::vector<SparseLinearClassifier*> classifiers;
    Vector w; ///< classifier weights
    Vector P; ///< classifier selection probabilities
    Vector output;
    real alpha;
    LinearClassifierMixture(int n_inputs_, int n_classes_,
                            int n_classifiers);
    virtual ~LinearClassifierMixture();
    int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector& Output(const Vector& x);
    virtual void Observe(const Vector& x, int label);
    void Show();
    void setStepSize(real step_size);
};

class HashedLinearClassifierMixture : public LinearClassifierMixture
{
protected:
    unsigned long secret;
public:
    
    HashedLinearClassifierMixture(int n_inputs_, int n_classes_, int n_classifiers); 
    virtual ~HashedLinearClassifierMixture()
    {
    }
    virtual void Observe(const Vector& x, int label);
};


template <class T>
class HashedClassifierMixture 
{
protected:
    const int n_inputs;
    const int n_classes;
    std::vector<T*> classifiers;
    unsigned long secret;
public:
    Vector output;
    HashedClassifierMixture(int n_inputs_,
                            int n_classes_,
                            std::vector<T*>& classifiers_)
        : n_inputs(n_inputs_),
          n_classes(n_classes_),
          classifiers(classifiers_),
          output(n_classes)
    {
        secret = true_random_bits(false);
    }
    ~HashedClassifierMixture()
    {
    }

    void Observe(Vector& x, int label)
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
            if (hash & 1) {
                //printf ("1");
                classifiers[i]->Observe(x, label);
            } else {
                //printf("0");
            }
            hash = hash >> 1;
        }
    }
    int Classify(Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector& Output(Vector& x)
    {
        output.Clear();
        real w = 1.0 / (real) classifiers.size();
        for (int i=0; i<(int) classifiers.size(); ++i) {
            output += classifiers[i]->Output(x);
        }
        output *= w;
        return output;
    }
};




#endif
