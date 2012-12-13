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

#include "BayesianMultivariate.h"
#include "BayesianMultivariateRegression.h"
#include "Random.h"
#include "Vector.h"

BayesianMultivariate::BayesianMultivariate(int n_actions_, 
										   int n_input_dim_,
										   int n_output_dim_,
										   real gamma_, 
										   real epsilon_,
										   RBFBasisSet* rbf_,
										   std::vector<BayesianMultivariateRegression*> regression_t_,
										   std::vector<BayesianMultivariateRegression*> regression_r_,
										   LinearModel<Vector,int>* lm_,
										   RepresentativeStateModel<LinearModel<Vector, int>, Vector, int>* RSM_,
										   real baseline_)
					: n_actions(n_actions_),
					  n_input_dim(n_input_dim_),
					  n_output_dim(n_output_dim_),
					  gamma(gamma_),
					  epsilon(epsilon_),
					  rbf(rbf_),
					  regression_t(regression_t_),
					  regression_r(regression_r_),
					  lm(lm_),
					  RSM(RSM_),
					  baseline(baseline_),
					  use_geometric_schedule(false)
{
	assert(n_grids > 0);
	assert( gamma >= 0.0 && gamma <= 1.0);
//  RepresentativeStateModel<LinearModel<Vector, int>, Vector, int> RR(gamma, *lm, n_states, n_actions);
//	real (RepresentativeStateModel<LinearModel<Vector, int>, Vector, int>::* getValue_)(const Vector&, const int&) = &RepresentativeStateModel<LinearModel<Vector, int>, Vector, int>::getValue; 
//	ContinuousStateEpsilonGreedy r(RR.getValue_, n_states, n_actions, epsilon);
//	Reset();
}

BayesianMultivariate::~BayesianMultivariate()
{	
}
void BayesianMultivariate::Reset()
{
	for(int i = 0; i<n_actions; ++i) {
		regression_t[i]->Reset();
		regression_r[i]->Reset();
	}
	lm->Reset();
	RSM->Reset();
}
void BayesianMultivariate::Observe(Vector state, int action, real reward, Vector next_state)
{	
	Vector phi;
	Vector r(reward);
	if(rbf != NULL) {
		rbf->Evaluate(state);
		phi = rbf->F();
//		phi.Resize(rbf->size() + 1);
//		phi[rbf->size()] = 1.0;
//		phi.print(stdout);
	}
	else {
		phi = state;
	}
	phi.Resize(n_input_dim);
	phi[n_input_dim - 1] = 1.0;
	
//	printf("Next State\n");
//	next_state.print(stdout);
//	printf("Prediction\n");
//	Vector out = regression_t[action]->generate()*phi;
//	out.print(stdout);
//	printf("Reward = %f\n",reward);
//	Vector a = regression_r[action]->generate()*phi;
//	printf("Predicted reward = %f\n",a[0]);
//	
	regression_t[action]->AddElement(next_state, phi);
	regression_r[action]->AddElement(r,phi);
	
//	printf("Next State\n");
//	next_state.print(stdout);
//	printf("Prediction\n");
//	Vector out = regression_t[action]->generate()*phi;
//	out.print(stdout);
//	printf("Reward = %f\n",reward);
//	Vector a = regression_r[action]->generate()*phi;
//	printf("Predicted reward = %f\n",a[0]);
	///////////
	//if (urandom() < (1 - gamma)) {
//		RSM->AddState(next_state);
//	}
	////////////
}
int BayesianMultivariate::Act(Vector state)
{
	Vector Q(n_actions);
	
	real threshold = epsilon;
	if(use_geometric_schedule) {
		threshold = epsilon / (1 + sqrt(beta));
		beta += alpha;
	}
	if(urandom() < threshold) {
		int action =  (int) floor(urandom(0.0, (real) Q.Size()));
		return action;
	}
	
	for (int i=0; i<Q.Size(); ++i) {
		Q(i) = RSM->getValue(state, i);
	}
	std::vector<int> max_elem = ArgMaxs(Q);
	return max_elem[urandom(0,max_elem.size())];	
}
real BayesianMultivariate::getValue(Vector state, int action)
{
	return RSM->getValue(state,action);
}
real BayesianMultivariate::getValue(Vector state)
{
	return RSM->getValue(state);
}
void BayesianMultivariate::Update()
{
    Matrix MeanR(1,n_input_dim);
	Matrix MeanT(n_output_dim, n_input_dim);
	Matrix CovarianceR(n_input_dim, n_input_dim);
	Matrix CovarianceT(n_input_dim, n_input_dim);
	
	for( int i=0; i<n_actions; ++i) {
		regression_t[i]->generate(MeanT, CovarianceT);
		lm->SetStatePredictionMean(MeanT,i);
		lm->SetStatePredictionVar(CovarianceT,i);
		
		regression_r[i]->generate(MeanR, CovarianceR);
		lm->SetRewardPredictionMean(MeanR.getRow(0),i);
		lm->SetRewardPredictionVar(CovarianceR,i);
	}
	RSM->Update(*lm);
}
void BayesianMultivariate::setGeometricSchedule(real alpha_, real beta_)
{
	alpha = alpha_;
	beta = beta_;
	use_geometric_schedule = true;
}
