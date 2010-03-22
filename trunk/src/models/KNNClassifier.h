/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef KNN_Classifier_H
#define KNN_Classifier_H

#include "Vector.h"
#include "KDTree.h"
#include <list>
#include <vector>


class KNNClassifier
{
public:
    class DataSample
    {
    public:
        Vector x; ///< point
        int label; ///< label
        DataSample(Vector& x_, int label_)
            : x(x_), label(label_)        
        {
        }
    };
protected:
    int n_classes;
    int n_dim;
    int K;
    KDTree<DataSample>* kd_tree;
    std::list<DataSample> samples;
    void AddSample(DataSample sample);
public:	
    Vector output;
    KNNClassifier(int n_classes_, int n_dim_, int K_);
    ~KNNClassifier();
    int Classify(Vector& x)
    {
        return ArgMax(Output(x));
    }
    Vector& Output(Vector& x);
    void Observe(Vector& x, int label)
    {
        AddSample(DataSample(x, label));
    }
};




#endif
