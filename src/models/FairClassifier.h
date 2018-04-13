/* -*- Mode: c++;  -*- */
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#if 0

#ifndef FAIR_CLASSIFIER_H
#define FAIR_CLASSIFIER_H

/** Generic fair classifier class.

	Given input \f$x\f$, sensitive variable \f$z\f$, and label
	\f$y\f$, build an appropriate classifier.

	Here we maintain a different classifier for each \f$z\f$, which gives us \fP(y, x, z) = P(y | x, z) P(x | z) P(z)\f$. We then feed this into the blind_classifier model.

	The number of classes may not be the same as the number of decisions. For example, when class may relate to a complicated merit variable, while the decisions may be simple 'accept/reject'.
 */
class FairClassifier
{
protected:
	std::vector<MultivariateGaussianClassifier*> classifier; ///< store \f$P(y | x, z)\f$.
	std::vector<MultivariateNormalUnknownMeanPrecision*> density; ///< store the distribution \f$P(x | z)\f$.
	LinearClassifier* blind_classifier; ///< Store \f$P(y | x)\f$.
	int n_inputs; ///< number of input dimensions
	int n_classes; ///< number of distinct classes
	int n_sensitive; ///< number of sensitive values
	int n_actions; ///< number of classification actions
	Dirichlet dirichlet; ///< stores the values of sensitive variables
	Matrix utility; ///< the utility of each decision-class pair
public:
	/// Create a fair classifier with n values of the sensitive variable
	FairClassifier(int n_inputs_, int n_classes_, int n_sensitive_, const Matrix& utility_);
    /// Destructor, where the classifiers are de-allocated
    virtual ~FairClassifier();
	/// Observe a new triplet of x, label and z
	virtual real Observe(const Vector& x,
						 const int label,
						 const int z);
    /// Classify a specific input according to our decision rule
    virtual int Classify(const Vector& x);
	/// Measure the utility of the classifiier for a given point
	virtual real Utility(const Vector& x, const int& z, const int& y);
	/// Measure the discrimination of a given policy
	virtual real Discrimination();
		
};

#endif

#endif
