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
										   FittedValueIteration<Vector,int>* FVI_,
										   FittedLSTD<Vector,int>* FLSTD_,
										   FittedQValueIteration<Vector,int>* FQVI_,
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
					  FVI(FVI_),
					  FLSTD(FLSTD_),
					  FQVI(FQVI_),
					  baseline(baseline_),
					  use_geometric_schedule(false)
{
	assert( gamma >= 0.0 && gamma <= 1.0);
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
	if(lm!=NULL) {
		lm->Reset();
		RSM->Reset();
	}
	else if( FVI!=NULL) {
		FVI->Reset();
	}
	else if(FLSTD != NULL) {
		return FLSTD->Reset();
	}
	else if(FQVI != NULL) {
		return FQVI->Reset();
	}
}
void BayesianMultivariate::Observe(Vector state, int action, real reward, Vector next_state)
{	
	Vector phi;
	Vector r(reward);
	if(rbf != NULL) {
		rbf->Evaluate(state);
		phi = rbf->F();
	}
	else {
		phi = state;
	}
	phi.Resize(n_input_dim);
	phi[n_input_dim - 1] = 1.0;

	regression_t[action]->AddElement(next_state, phi);
	regression_r[action]->AddElement(r,phi);

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
	real max = 0.0;
	if(RSM != NULL) {
		max = RSM->getValue(state,0);
	}
	else if(FVI != NULL) {
		max = FVI->getValue(state,0);
	}
	else if(FLSTD != NULL) {
		max = FLSTD->getValue(state,0);
	}
	else if(FQVI != NULL) {
		max = FQVI->getValue(state,0);
	}
	Q(0) = max;
	int select = 0;
	for(int i=1; i<Q.Size(); ++i) {
		if(RSM != NULL) {
			Q(i) = RSM->getValue(state,i);
		}
		else if(FVI != NULL) {
			Q(i) = FVI->getValue(state,i);
		}
		else if(FLSTD != NULL) {
			Q(i) = FLSTD->getValue(state,i);
		}
		else if(FQVI != NULL) {
			Q(i) = FQVI->getValue(state,i);
		}
		if(Q(i) > max) {
			max = Q(i);
			select = i;
		}
	}
//	Q.print(stdout);
//	for (int i=0; i<Q.Size(); ++i) {
//		Q(i) = RSM->getValue(state, i);
//	}
	
	//	std::vector<int> max_elem = ArgMaxs(Q);
	return select; //max_elem[urandom(0,max_elem.size()-1)];	
}

real BayesianMultivariate::getValue(Vector state, int action)
{
	if(RSM != NULL) {
		return RSM->getValue(state,action);
	}
	else if(FVI != NULL) {
		return FVI->getValue(state,action);
	}
	else if(FLSTD != NULL) {
		return FLSTD->getValue(state, action);
	}
	else if(FQVI != NULL) {
		return FQVI->getValue(state, action);
	}
	else {
		printf("Error\n");
		return -1;
	}
}

real BayesianMultivariate::getValue(Vector state)
{
	if(RSM != NULL) {
		return RSM->getValue(state);
	}
	else if(FVI != NULL) {
		return FVI->getValue(state);
	}
	else if(FLSTD != NULL) {
		return FLSTD->getValue(state);
	}
	else if(FQVI != NULL) {
		return FQVI->getValue(state);
	}
	else {
		printf("Error\n");
		return -1;
	}
}

void BayesianMultivariate::Update()
{
    Matrix MeanR(1,n_input_dim);
	Matrix MeanT(n_output_dim, n_input_dim);
	Matrix CovarianceR(n_input_dim, n_input_dim);
	Matrix CovarianceT(n_input_dim, n_input_dim);
	
	if(lm != NULL) {
		for( int i=0; i<n_actions; ++i) {
			regression_t[i]->generate(MeanT, CovarianceT);
			lm->SetStatePredictionMean(MeanT,i);
			lm->SetStatePredictionVar(CovarianceT,i);
			lm->GetStatePredictionVar(i).print(stdout);
			lm->SetStatePredictionS(regression_t[i]->getSxx(),i);

//			regression_r[i]->generate(MeanR, CovarianceR);
//			lm->SetRewardPredictionMean(MeanR.getRow(0),i);
//			lm->SetRewardPredictionVar(CovarianceR,i);
//			lm->SetRewardPredictionS(regression_r[i]->getSxx(),i);
		}
		RSM->Update(*lm);
	}
	else if(FVI != NULL) {
		FVI->Update(0.0001, 100);
	}
	else if(FLSTD != NULL) {
		FLSTD->Update(0.0001, 100);
	}
	else if(FQVI != NULL) {
		FQVI->Update(0.0001, 100);
	}
}
void BayesianMultivariate::Predict()
{
  int n;
	Vector r;
	Vector phi;
	Vector next_state;
	Vector state;
	char buffer[100];
	
	for( int i=0; i<n_actions; ++i) {
		Matrix mean_s = regression_t[i]->generate();
		Matrix mean_r = regression_r[i]->generate();
		n = sprintf(buffer, "Predicted_Output_samples_action_%d", i);
		FILE *output  = fopen(buffer,"w");
		n = sprintf(buffer, "Predicted_Reward_samples_action_%d", i);
		FILE *rew = fopen(buffer,"w");
		if(output!=NULL && rew !=NULL) {
			for( int s=0; s< RSM->getNSamples(); ++s) {
				state = RSM->getSample(s);
				if(rbf != NULL) {
					rbf->Evaluate(state);
					phi = rbf->F();
				}
				else {
					phi = state;
				}
				phi.Resize(n_input_dim);
				phi[n_input_dim - 1] = 1.0;
				next_state = mean_s*phi;
				next_state.print(output);
				r = mean_r*phi;
				r.print(rew);
			}
		}
		fclose(output);
		fclose(rew);
	}
}

void BayesianMultivariate::Predict(std::vector<Vector> samples)
{
	int n;
	Vector r;
	Vector phi;
	Vector next_state;
	Vector state;
	char buffer[100];
	for( int i=0; i<n_actions; ++i) {
		Matrix mean_s = regression_t[i]->generate();
		Matrix mean_r = regression_r[i]->generate();
		n = sprintf(buffer, "Predicted_Output_samples_action_%d", i);
		FILE *output  = fopen(buffer,"w");
		n = sprintf(buffer, "Predicted_Reward_samples_action_%d", i);
		FILE *rew = fopen(buffer,"w");
		if(output!=NULL && rew !=NULL) {
			for( uint s=0; s< samples.size(); ++s) {
				state = samples[s];
				if(rbf != NULL) {
					rbf->Evaluate(state);
					phi = rbf->F();
				}
				else {
					phi = state;
				}
				phi.Resize(n_input_dim);
				phi[n_input_dim - 1] = 1.0;
				next_state = mean_s*phi;
				next_state.print(output);
				r = mean_r*phi;
				r.print(rew);
			}
		}
		fclose(output);
		fclose(rew);
	}
	
	n = sprintf(buffer, "Predicted value function");
	FILE *value = fopen(buffer, "w");
	if(FVI != NULL) {
		if(value!=NULL) {
			Vector V((int)samples.size());
			for(uint s = 0; s < samples.size(); ++s) {
				V[s] = FVI->getValue(samples[s]);
			}
			V.print(value);
		}
	}
	else if(FLSTD != NULL) {
		if(value!=NULL) {
			Vector V((int)samples.size());
			for(uint s = 0; s < samples.size(); ++s) {
				V[s] = FLSTD->getValue(samples[s]);
			}
			V.print(value);
		}
	}
	else if(FQVI != NULL) {
		if(value!=NULL) {
			Vector V((int)samples.size());
			for(uint s = 0; s < samples.size(); ++s) {
				V[s] = FQVI->getValue(samples[s]);
			}
			V.print(value);
		}
	}
	fclose(value);
	
}

void BayesianMultivariate::setGeometricSchedule(real alpha_, real beta_)
{
	alpha = alpha_;
	beta = beta_;
	use_geometric_schedule = true;
}
