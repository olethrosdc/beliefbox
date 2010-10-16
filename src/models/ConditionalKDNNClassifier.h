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

#ifndef CONDITIONAL_KDNN_CLASSIFIER_H
#define CONDITIONAL_KDNN_CLASSIFIER_H


#include <vector>
#include "real.h"
#include "Vector.h"
#include "Ring.h"
#include "ContextTreeKDTree.h"
#include "NormalDistribution.h"
#include "KNNClassifier.h"
#include "Dirichlet.h"

/** Context tree non-parametric conditional density estimation on \f$R^n \times R^m\f$.

    This is a variant of the ContextTree class. The main difference is
    that it defines a conditional distribution \f$P(\cdot \mid x)\f$
    with \f$x \in R^n\f$ rather than the set of all strings
    \f$X^*\f$. However, \f$P\f$ remains a distribution on \f$X\f$.

	The model can also be used for estimating conditional
	distributions directly, and thus, for performing classification,
    directly.
 */
class ConditionalKDNNClassifier
{
public:
    // public classes
    struct Node
    {
		ConditionalKDNNClassifier& tree; ///< the tree
        Vector lower_bound_x; ///< looks at x > lower_bound
        Vector upper_bound_x; ///< looks at x < upper_bound
        real mid_point; ///< how to split
		int splitting_dimension; ///< dimension on which to do the split.
        const int depth; ///< depth of the node
        real log_prior_normal;
        real prior_normal;
		KNNClassifier* local_probability; ///< local probability estimator
        Node* prev; ///< previous node
        std::vector<Node*> next; ///< pointers to next nodes
        real w; ///< backoff weight
        real log_w; ///< log of w
		DirichletDistribution prior;
        Node(ConditionalKDNNClassifier& tree_,
			 Vector& lower_bound_x_,
			 Vector& upper_bound_x_);
        Node(Node* prev_, 
			 Vector& lower_bound_x, Vector& upper_bound_x);
        ~Node();
        real Observe(const Vector& x, const int y, real probability);
        real pdf(const Vector& x, const int y, real probability);
        Vector Output(const Vector& x, const Vector& P_y);
        void Show();
        int NChildren();    
        int S;
    };
    
    // public methods
    ConditionalKDNNClassifier(int n_branches_, int max_depth_,
                                    Vector& lower_bound_x, Vector& upper_bound_x,
                                    int n_classes_);
    ~ConditionalKDNNClassifier();
    real Observe(const Vector& x, const int y);
    int Classify(const Vector& x) 
    {
        return ArgMax(Output(x));
    }
    Vector& Output(const Vector& x);
    real pdf(const Vector& x, const int y);
    void Show();
    int NChildren();
protected: 
    int n_branches;
    int max_depth;
    int n_classes;
    Node* root;
    static real FUDGE;
public:
    Vector output;
};



#endif
