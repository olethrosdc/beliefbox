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

#ifndef KNN_CLASSIFIER_H
#define KNN_CLASSIFIER_H

#include "Vector.h"
#include "KDTree.h"
#include <list>
#include <vector>
#include "Classifier.h"


/** K Nearest neighbour classifier.
 */
class KNNClassifier : public  Classifier<Vector, int, Vector>
{
public:
    /// These are the samples stored in the KNN back-end.
    class DataSample
    {
    public:
        Vector x; ///< point
        ///int label; ///< label
        Vector Py; ///< class propabilities
        DataSample(const Vector& x_, const Vector& p)
            : x(x_), Py(p)        
        {
        }
    };
protected:
    int n_classes; ///< number of classes
    int n_dim; ///< number of input dimensions
    int K; ///< number of neighbours to use
    KDTree<DataSample> kd_tree; ///< storage backend
    std::list<DataSample> samples; ///< list of samples
    void AddSample(const DataSample sample); 
public:	
    Vector output; ///< temporary storage for the last output of the classifier
    KNNClassifier(const int n_dim_, const int n_classes_, const int K_);
    virtual ~KNNClassifier();
    virtual int Classify(const Vector& x)
    {
        return ArgMax(Output(x));
    }
    virtual Vector& Output(const Vector& x);
    virtual real Observe(const Vector& x, const int& label)
    {
		real p_y_x = Output(x)(label);
        Vector P(n_classes);
        P(label) = 1;
        AddSample(DataSample(x, P));
		return p_y_x;
    }
    virtual real Observe(const Vector& x, const Vector& p)
    {
        Vector p_y_x = Output(x);
        AddSample(DataSample(x, p));
		return Product(&p_y_x, &p);
    }
	void Show()
	{
		kd_tree.Show();
	}
};




#endif
