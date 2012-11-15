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
#include "Classifier.h"
#include "Vector.h"
#include <vector>

class LinearClassifierMixture : public Classifier<Vector, int, Vector>
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
                            int n_classifiers, int projection_size);
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
    
    HashedLinearClassifierMixture(int n_inputs_, int n_classes_, int n_classifiers, int projection_size); 
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
	Vector P;
	Vector w;
public:
    Vector output;
    HashedClassifierMixture(int n_inputs_,
                            int n_classes_,
                            std::vector<T*>& classifiers_)
        : n_inputs(n_inputs_),
          n_classes(n_classes_),
          classifiers(classifiers_),
		  P((int) classifiers.size()),
		  w((int) classifiers.size()),
          output(n_classes)		  
    {
		real w_prior = 1.0 / (real) w.Size();
		for (int i=0; i<w.Size(); ++i) {
			w(i) = w_prior;
		}
        secret = true_random_bits(false);
    }
    ~HashedClassifierMixture()
    {
    }

    void Observe(Vector& x, int label)
    {
        assert(x.Size() == n_inputs);
        assert(label >= 0 && label < n_classes);
    
		for (int i=0; i<(int) classifiers.size(); ++i) {
			P(i) = w(i) * (classifiers[i]->Output(x))(label);
		}
		w = P / P.Sum();
    
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
    void Show()
	{
#if 0
		for (int i=0; i<(int) classifiers.size(); ++i) {
			printf ("# %d\n", i);
			classifiers[i]->Show();
		}
#endif
		printf ("w: ");
		w.print(stdout);
	}
};




#endif
