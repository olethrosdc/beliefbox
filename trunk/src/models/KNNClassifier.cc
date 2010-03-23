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

#include "KNNClassifier.h"
#include "BasisSet.h"
#include "Distribution.h"
#include "Random.h"

/** Create a model
    
    \param n_classes_ set the number of classes
    \param n_dim_ the number of dimensions in the observation space
    \param K the number of neighbours
    
    As a side-effect, initialise the kd_tree list
 */
KNNClassifier::KNNClassifier(int n_classes_, int n_dim_, int K_)
    : n_classes(n_classes_), n_dim(n_dim_), K(K_),
      kd_tree(n_dim), output(n_classes)
{
}

/// Delete model, as well as the kd_tree
KNNClassifier::~KNNClassifier()
{
}

/** Add a sample to the model
    
    \param sample sample to add
    
    Add a sample if not too many samples have been added.  The sample
    is added to the list of samples.  It is also added to the tree,
    using x as a context.
    
 */
void KNNClassifier::AddSample(DataSample sample)
{
    samples.push_back(sample);
    kd_tree.AddVectorObject(sample.x, &samples.back());
}


/** Predict the next label.
    
    This is the simplest type of predictor.
    It just predicts the next vector mean, no distribution is used.
    
    \param x observables
*/
Vector& KNNClassifier::Output(Vector& x)
{
    //basis.Evaluate(x);
    assert(n_dim == x.Size());
    assert(n_classes == output.Size());
    for (int i=0; i<n_classes; ++i) {
        output[i] = 0;
    }

    real w = 1.0 / (real) K;

    OrderedFixedList<KDNode> node_list = kd_tree.FindKNearestNeighbours(x, K);
    std::list<std::pair<real, KDNode*> >::iterator it;
    int i=0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it, ++i) {
        KDNode* node = it->second;
        DataSample* sample = kd_tree.getObject(node);
        output[sample->label] += w;
    }
    return output;
}




