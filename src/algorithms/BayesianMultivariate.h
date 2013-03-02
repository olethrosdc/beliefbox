// -*- Mode: c++ -*-
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BAYESIAN_MULTIVARIATE_H
#define BAYESIAN_MULTIVARIATE_H

#include "DiscretePolicy.h"
#include "BasisSet.h"
#include "RepresentativeStateModel.h"
#include "BayesianMultivariateRegression.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "ExplorationPolicy.h"
#include <vector>
#include "LinearModel.h"
#include "FittedValueIteration.h"
#include "FittedQValueIteration.h"
#include "FittedLSTD.h"

class BayesianMultivariate 
{
protected:
	const int n_actions;	///< total number of actions
	const int n_input_dim;	///< Input dimension
	const int n_output_dim; ///< Output dimension
	real gamma;				///< discount factor
//	ContinuousStateEpsilonGreedy* exploration_policy; ///<exploration policy
	real epsilon;
	RBFBasisSet* rbf;
	std::vector<BayesianMultivariateRegression*> regression_t;
	std::vector<BayesianMultivariateRegression*> regression_r;
	LinearModel<Vector, int>* lm;
	RepresentativeStateModel<LinearModel<Vector, int>, Vector, int>* RSM;
	FittedValueIteration<Vector,int>* FVI;
	FittedLSTD<Vector,int>* FLSTD;
	FittedQValueIteration<Vector, int>* FQVI;
	real baseline;
	real beta;
	real alpha;
    bool use_geometric_schedule;
public:
	BayesianMultivariate(int n_actions_, 
						 int n_input_dim_,
						 int n_output_dim_,
						 real gamma_, 
						 real epsilon_,
						 RBFBasisSet* rbf_,
						 std::vector<BayesianMultivariateRegression*> regression_t_,
						 std::vector<BayesianMultivariateRegression*> regression_r_,
						 LinearModel<Vector,int>* lm_ = NULL,
						 RepresentativeStateModel<LinearModel<Vector, int>, Vector, int>* RSM_ = NULL,
						 FittedValueIteration<Vector, int>* FVI_ = NULL,
						 FittedLSTD<Vector,int>* FLSTD_ = NULL,
						 FittedQValueIteration<Vector, int>* FQVI_ = NULL,
						 real baseline_ = 0.0);
	virtual ~BayesianMultivariate();
	virtual void Reset();
	/// Full bayesian multivariate observation.
	virtual void Observe(Vector state, int action, real reward, Vector next_state);
	///
	virtual int Act(Vector next_state);
	///
	virtual real getValue(Vector state, int action);
	///
	real getValue(Vector state);
	///
	virtual void Update();
	///
	void Predict();
	///
	void Predict(std::vector<Vector> samples);
	///
	void setGeometricSchedule(real alpha_, real beta_);
    ///
};

#endif