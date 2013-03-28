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

#include "CoverBayesianMultivariate.h"
#include "CoverTree.h"
#include "Random.h"
#include "Vector.h"

CoverBayesianMultivariate::CoverBayesianMultivariate(int n_actions_, 
										   real gamma_, 
										   real epsilon_,
										   std::vector<CoverTree*> cover_,
										   CoverFittedValueIteration<Vector, int>* FVI_,
										   CoverFittedLSTD<Vector,int>* FLSTD_,
										   CoverFittedLSTDQ<Vector, int>* FLSTDQ_,
										   real baseline_)
		: n_actions(n_actions_),
		gamma(gamma_),
		epsilon(epsilon_),
		cover(cover_),
		FVI(FVI_),
		FLSTD(FLSTD_),
		FLSTDQ(FLSTDQ_),
		baseline(baseline_),
		use_geometric_schedule(false) 
{
	assert( gamma >= 0.0 && gamma <= 1.0);
	assert( epsilon >= 0 );
}

CoverBayesianMultivariate::~CoverBayesianMultivariate()
{	
}

void CoverBayesianMultivariate::Reset()
{
	for(int i = 0; i<n_actions; ++i) {
		cover[i]->Reset();
	}
	if( FVI != NULL) {
		FVI->Reset();
	}
	else if(FLSTD != NULL) {
		return FLSTD->Reset();
	}
	else if(FLSTDQ != NULL) {
		return FLSTDQ->Reset();
	}
	else {
		printf("You haven't declare the FittedValueIteration algorithm\n");
	}
	//if(lm!=NULL) {
//		lm->Reset();
//		RSM->Reset();
//	}
//	else if( FVI!=NULL) {
//		FVI->Reset();
//	}
//	else if(FLSTD != NULL) {
//		return FLSTD->Reset();
//	}
//	else if(FQVI != NULL) {
//		return FQVI->Reset();
//	}
}
void CoverBayesianMultivariate::Observe(Vector state, int action, real reward, Vector next_state)
{	
	//CoverTree::Node* node = 
	cover[action]->Insert(state, next_state, reward, false);

//	CoverTree::Node* node = cover[action]->Insert(state);
}
int CoverBayesianMultivariate::Act(Vector state)
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
	
	if(FVI != NULL) {
		max = FVI->getValue(state, 0);
	}
	else if(FLSTD != NULL) {
		max = FLSTD->getValue(state,0);
	}
	else if(FLSTDQ != NULL) {
		max = FLSTDQ->getValue(state,0);
	}
	
	Q(0) = max;
	int select = 0;
	for(int i=1; i<Q.Size(); ++i) {
		if(FVI != NULL) {
			Q(i) = FVI->getValue(state,i);
		}
		else if(FLSTD != NULL) {
			Q(i) = FLSTD->getValue(state,i);
		}
		else if(FLSTDQ != NULL) {
			Q(i) = FLSTDQ->getValue(state,i);
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
real CoverBayesianMultivariate::getValue(Vector state, int action)
{
	if( FVI != NULL) {
		return FVI->getValue(state, action);
	}
	else if(FLSTD != NULL) {
		return FLSTD->getValue(state, action);
	}
	else if(FLSTDQ != NULL) {
		return FLSTDQ->getValue(state, action);
	}
	else {
		printf("Error\n");
		return -1;
	}
}
real CoverBayesianMultivariate::getValue(Vector state)
{
	if( FVI != NULL) {
		return FVI->getValue(state);
	}
	else if(FLSTD != NULL) {
		return FLSTD->getValue(state);
	}	
	else if(FLSTDQ != NULL) {
		return FLSTDQ->getValue(state);
	}
	else {
		printf("Error\n");
		return -1;
	}
}
void CoverBayesianMultivariate::Update()
{
	for( int i=0; i<n_actions; ++i) {
		cover[i]->SamplingTree();
	} 
	if(FVI != NULL) {
		FVI->Update(0.0001, 10);
	}
	else if(FLSTD != NULL) {
		FLSTD->Update(0.0001, 20);
	}
	else if(FLSTDQ != NULL) {
		FLSTDQ->Update(0.0001, 10);
	}
	else {
		printf("Problem on CoverBayesianMultivariate\n");
	}
}
//void CoverBayesianMultivariate::Predict()
//{
//	int n;
//	Vector r;
//	Vector phi;
//	Vector next_state;
//	Vector state;
//	char buffer[100];
//	
//	for( int i=0; i<n_actions; ++i) {
//		Matrix mean_s = regression_t[i]->generate();
//		Matrix mean_r = regression_r[i]->generate();
//		n = sprintf(buffer, "Predicted_Output_samples_action_%d", i);
//		FILE *output  = fopen(buffer,"w");
//		n = sprintf(buffer, "Predicted_Reward_samples_action_%d", i);
//		FILE *rew = fopen(buffer,"w");
//		if(output!=NULL && rew !=NULL) {
//			for( int s=0; s< RSM->getNSamples(); ++s) {
//				state = RSM->getSample(s);
//				if(rbf != NULL) {
//					rbf->Evaluate(state);
//					phi = rbf->F();
//				}
//				else {
//					phi = state;
//				}
//				phi.Resize(n_input_dim);
//				phi[n_input_dim - 1] = 1.0;
//				next_state = mean_s*phi;
//				next_state.print(output);
//				r = mean_r*phi;
//				r.print(rew);
//			}
//		}
//		fclose(output);
//		fclose(rew);
//	}
//}

void CoverBayesianMultivariate::Predict(std::vector<Vector> samples)
{
	int n;
	Vector phi;
	Vector next_state;
	Vector state;
	char buffer[100];
	
	for( int i=0; i<n_actions; ++i) {
		n = sprintf(buffer, "Predicted_Output_samples_action_%d", i);
		cover[i]->SamplingTree();
		FILE *output  = fopen(buffer,"w");
		if(output!=NULL) {
			for( uint s=0; s< samples.size(); ++s) {
				next_state = cover[i]->GenerateState(samples[s]);
				next_state.print(output);
			}
		}
		fclose(output);
	}
	
	n = sprintf(buffer, "Predicted value function");
	FILE *value = fopen(buffer, "w");
	if(FLSTD != NULL) {
		if(value!=NULL) {
			Vector V((int)samples.size());
			for(uint s = 0; s < samples.size(); ++s) {
				V[s] = FLSTD->getValue(samples[s]);
			}
			V.print(value);
		}
	}
	fclose(value);
	
}

void CoverBayesianMultivariate::setGeometricSchedule(real alpha_, real beta_)
{
	alpha = alpha_;
	beta = beta_;
	use_geometric_schedule = true;
}
