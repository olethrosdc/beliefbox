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

#include "ConditionalKDGaussianClassifier.h"
#include "BetaDistribution.h"
#include <cmath>

#include "Random.h"

ConditionalKDGaussianClassifier::Node::Node(ConditionalKDGaussianClassifier& tree_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(tree_),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(0),
      prev(NULL),
      next(tree.n_branches),
      w(1),
	  log_w(0),
      S(0)
{
    assert(lower_bound_x < upper_bound_x);
    //assert(lower_bound_y < upper_bound_y);
	splitting_dimension = -1 ;//ArgMax(upper_bound_x - lower_bound_x);
#ifdef RANDOM_SPLITS
	real phi = 0.1 + 0.8 * urandom();
	mid_point = phi * upper_bound_x[splitting_dimension] + 
		(1 - phi) * lower_bound_x[splitting_dimension];
#else
	mid_point = 0;//(upper_bound_x[splitting_dimension] + lower_bound_x[splitting_dimension]) / 2.0;
#endif


	local_probability = new MultivariateGaussianClassifier(lower_bound_x.Size(), tree.n_classes);
}

/// Make a node for K symbols at nominal depth d
ConditionalKDGaussianClassifier::Node::Node(ConditionalKDGaussianClassifier::Node* prev_,
									 Vector& lower_bound_x_,
									 Vector& upper_bound_x_)
    : tree(prev_->tree),
	  lower_bound_x(lower_bound_x_),
      upper_bound_x(upper_bound_x_),
      depth(prev_->depth + 1),
      prev(prev_),
      next(tree.n_branches),
      log_w(- depth * log(2)),
      S(0)
{
    assert(lower_bound_x < upper_bound_x);
	splitting_dimension = -1;//ArgMax(upper_bound_x - lower_bound_x);    
#ifdef RANDOM_SPLITS
	real phi = 0.1 + 0.8 * urandom();
	mid_point = phi * upper_bound_x[splitting_dimension] + 
		(1 - phi) * lower_bound_x[splitting_dimension];
#else
	mid_point = 0;//(upper_bound_x[splitting_dimension] + lower_bound_x[splitting_dimension]) / 2.0;
#endif

    w = exp(log_w);
    for (int i=0; i<tree.n_branches; ++i) {
        next[i] = NULL;
    }

	local_probability = new MultivariateGaussianClassifier(lower_bound_x.Size(),
                                                           tree.n_classes);
}

/// make sure to kill all
ConditionalKDGaussianClassifier::Node::~Node()
{
	delete local_probability;
    for (int i=0; i<tree.n_branches; ++i) {
        delete next[i];
    }
}

/** Observe new data, adapt parameters.
    
    Each node corresponds to an interval C.
	We only want to model the conditional $\Pr(y \mid x \in C)\f$.
*/
real ConditionalKDGaussianClassifier::Node::Observe(const Vector& x, const int y, real probability)
{
	real total_probability;

    // the local distribution
    real P_local = 0.001 / (real) tree.n_classes + 0.999 * local_probability->Observe(x, y);

	// Mixture with the previous ones
    w = exp(log_w); 
    assert (w >= 0 && w <= 1);
	total_probability = P_local * w + (1 - w) * probability;

    // adapt parameters
    log_w += log(P_local) - log(total_probability);
    w = exp(log_w); 
    assert (w >= 0 && w <= 1);


    real threshold = depth; //sqrt((real) depth);
    S++;

    DirichletDistribution& prior = local_probability->prior;
    int most_observed_class = ArgMax(prior.GetMean());
    //real n_most_observed_class = prior.GetMean()(most_observed_class);
    int second_most_observed_class = -1;
    real n_second_most_observed_class = 0;
    int number_of_observed_classes = 0;
    Vector A=prior.GetParameters();
    for (int i=0; i<tree.n_classes; ++i) {
        if (A(i) > 1) {
            number_of_observed_classes++;
            if (i != most_observed_class 
                && A(i) > n_second_most_observed_class) {
                second_most_observed_class = i;
                n_second_most_observed_class = A(i);
            }
        }
    }

#if 0
    std::cout << "observed classes: " << number_of_observed_classes
              << ", most observed: " << most_observed_class
              << ", second most observed: " << second_most_observed_class
              << std::endl;

    std::cout << depth << ": P(y|h_k)=" << P_local
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << std::endl;
#endif
	// Do a forward mixture if there is another node available.
    if ((tree.max_depth<0 || depth < tree.max_depth) && S >  threshold && number_of_observed_classes > 1) {
        if (splitting_dimension < 0) {
            Vector mean = local_probability->getClassMean(most_observed_class);
            Vector other_mean = local_probability->getClassMean(second_most_observed_class);
            splitting_dimension = ArgMax(abs(mean - other_mean));
            mid_point = 0.5 * (mean(splitting_dimension) + other_mean(splitting_dimension));
            if (mid_point > upper_bound_x(splitting_dimension)
                || mid_point < lower_bound_x(splitting_dimension)) {
                std::cerr << "Warning: mid point exceeds bounds\n";
                mid_point = 0.5*(upper_bound_x(splitting_dimension) + lower_bound_x(splitting_dimension));
            }
                
#if 0
            std::cout << "# d: " << depth << ", i: " 
                      << splitting_dimension
                      << " x_i: " << mid_point 
                      << "[" << mean(splitting_dimension) 
                      << ", " << other_mean(splitting_dimension) << "]"
                      << std::endl;
#endif
        }
        // Which interval is the x lying at
        int interval = 0;
        if (splitting_dimension >= 0 && x[splitting_dimension] < mid_point) {
            interval = 0;
        } else {
            interval = 1;
        }
        if (!next[interval]) {

            if (interval == 0) {
				Vector new_bound_x = upper_bound_x;
				new_bound_x(splitting_dimension) = mid_point;
                next[interval] = new Node(this, lower_bound_x, new_bound_x);
            } else {
				Vector new_bound_x = lower_bound_x;
				new_bound_x(splitting_dimension) = mid_point;
                next[interval] = new Node(this, new_bound_x, upper_bound_x);
            }
        }
		total_probability = next[interval]->Observe(x, y, total_probability);
    }

    return total_probability;
}

/** Recursively obtain \f$P(y | x)\f$.
*/
real ConditionalKDGaussianClassifier::Node::pdf(const Vector& x, const int y, real probability)
{
	real total_probability;

    // the local distribution
    real P_local =  0.001 / (real) tree.n_classes + 0.999 * local_probability->Output(x)(y);

	// Mix the current one with all previous ones
    w =  exp(log_w);
	total_probability = P_local * w + (1 - w) * probability;

#if 0
    std::cout << depth << ": P(y|h_k)=" << P_local
              << ", P(h_k|B_k)=" << w 
              << ", P(y|B_{k-1})="<< probability
              << ", P(y|B_k)=" << total_probability
              << ", P(N)= " << prior_normal
              << std::endl;
#endif
    // Which interval is the x lying at
    if (splitting_dimension >= 0) {
        int k;
        if ( x[splitting_dimension] < mid_point) {
            k = 0;
        } else {
            k = 1;
        }    
        // Do one more mixing step if required
        if (next[k]) {
			total_probability = next[k]->pdf(x, y, total_probability);
        }
    }
    return total_probability;
}

/** Recursively obtain \f$P(y | x)\f$.
*/
Vector ConditionalKDGaussianClassifier::Node::Output(const Vector& x, const Vector& P_y)
{
	Vector total_probability(tree.n_classes);

    // the local distribution
    Vector P_local =   local_probability->Output(x) * 0.999 + 0.001 / (real) tree.n_classes;

	// Mix the current one with all previous ones
    w =  exp(log_w);
	total_probability = P_local * w + P_y * (1 - w);
#if 0
    std::cout << depth << "= ";
    P_local.print(stdout);
    std::cout << depth << "< ";
    P_y.print(stdout);
    std::cout << depth 
              << ", P(h_k|B_k)=" << w 
              << std::endl;
#endif
    // Which interval is the x lying at
    if (splitting_dimension >= 0) {
        int k;
        if ( x[splitting_dimension] < mid_point) {
            k = 0;
        } else {
            k = 1;
        }    
        // Do one more mixing step if required
        if (next[k]) {
			total_probability = next[k]->Output(x, total_probability);
        }
    }
    return total_probability;
}



void ConditionalKDGaussianClassifier::Node::Show()
{
#if 1
	printf("%d %f #w\n", depth, w);
	for (int k=0; k<tree.n_branches; ++k) {
		if (next[k]) {
			next[k]->Show();
		}
	}
#endif
	return;
}

int ConditionalKDGaussianClassifier::Node::NChildren()
{
    int my_children = 0;
    for (int k=0; k<tree.n_branches; ++k) {
        if (next[k]) {
            my_children++;
            my_children += next[k]->NChildren();
        }
    }
    return my_children;
}

ConditionalKDGaussianClassifier::ConditionalKDGaussianClassifier(int n_branches_,
                                                                 int max_depth_,
                                                                 Vector& lower_bound_x,
                                                                 Vector& upper_bound_x,
                                                                 int n_classes_)
    : n_branches(n_branches_),
      max_depth(max_depth_),
      n_classes(n_classes_),
      output(n_classes)
{
    root = new Node(*this, lower_bound_x, upper_bound_x);
}

ConditionalKDGaussianClassifier::~ConditionalKDGaussianClassifier()
{
    delete root;
}

/** Obtain \f$\xi_t(y \mid x)\f$ and calculate \f$\xi_{t+1}(w) = \xi_t(w \mid x, y)\f$.
 */
real ConditionalKDGaussianClassifier::Observe(const Vector& x, const int y)
{
    return root->Observe(x, y, 1);
}

/** Obtain \f$\xi_t(y \mid x)\f$ */
real ConditionalKDGaussianClassifier::pdf(const Vector& x, const int y)
{
    return root->pdf(x, y, 1);
}

/** Obtain \f$\xi_t(y \mid x)\f$ */
Vector& ConditionalKDGaussianClassifier::Output(const Vector& x)
{
    Vector tmp(n_classes);
    output = root->Output(x, tmp);
    return output;
}

void ConditionalKDGaussianClassifier::Show()
{
    root->Show();
    std::cout << "Total contexts: " << NChildren() << std::endl;
}

int ConditionalKDGaussianClassifier::NChildren()
{
    return root->NChildren();
}
