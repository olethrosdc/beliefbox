
#include "FairClassifier.h"

FairClassifier::FairClassifier(int n_inputs_,
							   int n_classes_,
							   int n_sensitive_,
							   const Matrix& utility_)
	: n_inputs(n_classes_),
	  n_classes(n_classes_),
	  n_sensitive(n_sensitive_),
	  dirichlet(n_sensitive),
	  utility(utility_),
	  n_actions(utility.Rows());
{
	classifier.resize(n_sensitive);
	for (int i=0; i<n_sensitive; i++) {
		classifier[i] = new MultivariateGaussianClassifier(n_inputs, n_classes);
		density[i] = new MultivariateNormalUnknownMeanPrecision(n_inputs);
	}
	blind_classifier = new MultivariateGaussianClassifier(n_inputs, n_classes);
	
}

FairClassifier::~FairClassifier()
{
	for (int i=0; i<n_sensitive; i++) {
		delete classifier[i];
		delete density[i];
	}
	delete blind_classifier;
}

real FairClassifier::Observe(const Vector& x,
									 const int label,
									 const int z)
{
	/// TO-DO
}


int FairClassifier::Classify(const Vector& x)
{
	for (int z=0; z<n_sensitive; z++) {
		
	}
}


































