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

#ifndef COVER_BAYESIAN_MULTIVARIATE_H
#define COVER_BAYESIAN_MULTIVARIATE_H

#include "DiscretePolicy.h"
#include "BasisSet.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "ExplorationPolicy.h"
#include <vector>
#include "LinearModel.h"
#include "CoverTree.h"
#include "CoverFittedValueIteration.h"
#include "CoverFittedLSTD.h"
#include "CoverFittedLSTDQ.h"

class CoverBayesianMultivariate 
{
	protected:
		const int n_actions;	///< total number of actions
		real gamma;				///< discount factor
		//	ContinuousStateEpsilonGreedy* exploration_policy; ///<exploration policy
		real epsilon;
		std::vector<CoverTree*> cover;
		CoverFittedValueIteration<Vector,int>* FVI;
		CoverFittedLSTD<Vector,int>* FLSTD;
		CoverFittedLSTDQ<Vector, int>* FLSTDQ;
		real baseline;
		real beta;
		real alpha;
		bool use_geometric_schedule;
	public:
		CoverBayesianMultivariate(int n_actions_, 
							 real gamma_, 
							 real epsilon_,
							 std::vector<CoverTree*> cover_,
							 CoverFittedValueIteration<Vector, int>* FVI_ = NULL,
							 CoverFittedLSTD<Vector,int>* FLSTD_ = NULL,
							 CoverFittedLSTDQ<Vector, int>* FLSTDQ_ = NULL,
							 real baseline_ = 0.0);
		virtual ~CoverBayesianMultivariate();
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
		void Predict(std::vector<Vector> samples, bool LinearTest = false);
		///
		void setGeometricSchedule(real alpha_, real beta_);
		///
	};

#endif
