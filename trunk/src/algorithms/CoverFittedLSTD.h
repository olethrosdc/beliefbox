/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef COVER_FITTED_LSTD_VALUE_ITERATION_H
#define COVER_FITTED_LSTD_VALUE_ITERATION_H

#include "BayesianMultivariateRegression.h"
#include "BasisSet.h"
#include "Environment.h"
#include "Matrix.h"
#include "Random.h"
#include "Rollout.h"
#include "MersenneTwister.h"
#include "RandomPolicy.h"
#include "ContinuousPolicy.h"
#include "TilingSet.h"
#include "vector.h"
#include <cstring>

///*Fitted Value Iteration Algorithm*/
template <class S, class A>
class CoverFittedLSTD {
protected:
    real gamma;					///< discount factor
	int N;						///< Number of basis points
	int M;						///< Number of sampled next states
	int dim;					///< Dimension of the basis functions
	int n_actions;				///< Number of actions (Note: what to do for continuous?)
	Vector weights;				///< Model parameters
	Matrix PHI;
	RBFBasisSet* RBFs;		
	real lambda;				///< Regularization factor
	Matrix pseudo_inv;
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<CoverTree*> cover;
	int grids;					///< Number of grids of basis function
	real scale;
	int tilings;				///< Number of tilings
	bool update_samples;		///< Take new random training samples 
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
public:
	CoverFittedLSTD(const real& gamma_, const int& N_, const int& M_, Environment<S,A>* environment_, std::vector<CoverTree*> cover_, const int& grids_  = 0, real scale_ = 1.0, const int& tilings_ = 0, bool update_samples_ = false):
	gamma(gamma_),
	N(N_),
	M(M_),
	environment(environment_),
	cover(cover_),
	grids(grids_),
	scale(scale_),
	tilings(tilings_),
	update_samples(update_samples_)
	{
		if(!strcmp(environment->Name(), "Bike")) {
			dim = 20;
		}
		else if(grids == 0) {
			RBFs = NULL;
			dim = 0;
			for( int i=0; i<n_actions; ++i) {
				dim += cover[i]->GetNumBasisNodes();
//				dim += cover[i]->GetNumSamplingNodes();
//				dim += cover[i]->GetNumNodes();
			}
		}
		else {
			if(tilings == 0) {
				dim	= pow(grids,(real)environment->getNStates()) + 1;	
				Vector S_L	= environment->StateLowerBound();
				Vector S_U	= environment->StateUpperBound();
				EvenGrid Discretisation(S_L,S_U,grids);
				RBFs = new RBFBasisSet(Discretisation,scale); 
			}
			else if(tilings > 0) {
				dim = 5*tilings*pow(grids,(real)environment->getNStates());
				printf("# %d Tilings (%d x %d) was created (Total basis functions =%d)\n",tilings, grids,grids,dim);
			}
		}
//		lambda = 0.01;
		lambda = 0.01;
		n_actions = environment->getNActions();
		
		sampleSelection();
		weights = Vector(dim);
	}
	const void Update(real threshold = 0.00001, int max_iter = -1) {
		
//		Reset();
		
		int a;
	    Matrix AA;
		Vector b(dim);
		Vector dif;
		Vector phi_;
		Vector phi;
		
		real r;
		bool final;
		real distance;
		Vector T(N);
		Vector pW(dim);
		weights = Vector::Vector(dim); // + 10;
		Vector V(n_actions);
		int n_iter = 0;
		do {
			AA = Matrix::Null(dim,dim);
			b.Clear();
			distance = 0.1;
			real steps = 0.0;
			for(int i=0; i<N; ++i)
			{	
//				printf("#%d\n",i);
				//Find best action - policy.
				for(a = 0; a < n_actions; a++) {
					environment->Reset();
					environment->setState(states[i]);
					
					final = environment->Act(a);
					r = environment->getReward();
					Vector next_state = cover[a]->GenerateState(states[i]);
					V[a] = getValue(next_state) + r;
//					printf("V[%d] => %f and reward = %f\n",a, V[a], r);
				}
		//		std::vector<int> v;
//				for(int aa = 0; aa < n_actions; aa++) {
//					if(V[aa] == Max(V)) {
//						v.push_back(aa);
//					}
//				}
//				int ss = urandom(0, (int)(v.size()));
//				a = v[ss];
				a = ArgMax(V);
				
				environment->Reset();
				environment->setState(states[i]);
				
				phi = BasisConstruction(states[i]);
//				phi.print(stdout);
				final = environment->Act(a);
				r = environment->getReward();
				
				if(final) {
					Vector next_state = cover[a]->GenerateState(states[i]);
					phi_ = BasisConstruction(next_state);
					dif = phi - (phi_*gamma);
				}
				else {
					dif = phi;
				}	
				steps++;
				Matrix res = OuterProduct(phi,dif);
				AA += res;					
				b  += phi*r;
			}
//			for(int i=0; i<rollout->getNRollouts(); ++i) {
//				for(int j=0; j<rollout->getNSamples(i); ++j) {
//					//				printf("#%d\n",i);
//					//Find best action - policy.
//					for(a = 0; a < n_actions; a++) {
//						environment->Reset();
//						environment->setState(rollout->getNextState(i,j));
//					
//						final = environment->Act(a);
//						r = environment->getReward();
////						Vector next_state = cover[a]->GenerateState(rollout->getNextState(i,j));
//						Vector next_state = environment->getState();
//						V[a] = getValue(next_state) + r;
//						//	printf("V[%d] => %f and reward = %f\n",a, V[a], r);
//					}
//				
//					a = ArgMax(V);
//				
//					environment->Reset();
//					environment->setState(rollout->getNextState(i,j));
//				
//					phi = BasisConstruction(rollout->getNextState(i,j));
//					//				phi.print(stdout);
//					final = environment->Act(a);
//					r = environment->getReward();
//				
//					if(final) {
//						Vector next_state = environment->getState();
////						Vector next_state = cover[a]->GenerateState(rollout->getNextState(i,j));
//						phi_ = BasisConstruction(next_state);
//						dif = phi - (phi_*gamma);
//					}
//					else {
//						dif = phi;
//					}	
//					steps++;
//					Matrix res = OuterProduct(phi,dif);
//					AA += res;					
//					b  += phi*r;
//				}
//			}
//			AA.print(stdout);
			weights = ((1.0/steps)*AA + (lambda*steps)*Matrix::Unity(dim,dim)).Inverse_LU()*(b*(1/steps));
			weights.print(stdout);
			distance = (pW - weights).L2Norm();
			pW = weights;
			if(max_iter>0) 
				max_iter--;
			n_iter++;
			
			if(update_samples==true) {
				sampleSelection();
			}
			printf("#ValueIteration (Threshold = %f) :: ComputeStateValues Exiting at d:%f, n:%d\n", threshold, distance, n_iter);
		}while(distance > threshold && max_iter != 0);
		printf("#ValueIteration::    ComputeStateValues Exiting at d:%f, n:%d\n", distance, n_iter);
	}
	
	real getValue(const S& state) 
	{
		Vector phi = BasisConstruction(state);	

		return Product(weights,phi);
	}
	
	real getValue(const S& state, const A& action) 
	{
		bool endsim = environment->getEndsim();
		Vector true_state = environment->getState();
		environment->Reset();
		environment->setState(state);
		environment->Act(action);
		real r = environment->getReward();
		
		environment->Reset();
		environment->setEndsim(endsim);
		environment->setState(true_state);
		
		real temp_v = 0.0;
		Vector next_state = cover[action]->GenerateState(state);
		temp_v = getValue(next_state);
		
		return (r + gamma*temp_v);
	}
	
	void sampleSelection() {
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
//		PHI = Matrix(dim,N);
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
//			Vector phi = BasisConstruction(state);			
//			PHI.setColumn(i,phi);
		}
//		pseudo_inv = (PHI*Transpose(PHI) + lambda*Matrix::Unity(dim,dim)).Inverse()*PHI;
	}
	
//	void sampleSelection() {
//		int n_rollouts	= 1000;
//		int horizon		= 20;
//		
//		RandomNumberGenerator* rng;
//		RandomNumberGenerator* environment_rng;
//		
//		MersenneTwisterRNG mersenne_twister_env;
//		environment_rng = (RandomNumberGenerator*) &mersenne_twister_env;
//		
//		MersenneTwisterRNG mersenne_twister;
//		rng = (RandomNumberGenerator*) &mersenne_twister;
//		
//		rng->manualSeed(123456789);
//		// Place holder for the policy
//        AbstractPolicy<Vector, int>* policy;
//		
//        // Start with a random policy!
//        policy = new RandomPolicy(environment->getNActions(), rng);
//		
//		environment->Reset();
//		Vector state_init = environment->getState();
//        rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> >(state_init, policy, environment, gamma, true);
//		
//		//		rollout->UniformSampling(n_rollouts, horizon);
//		//        rollout->StartingDistributionSampling(n_rollouts, horizon);
//		rollout->Sampling(n_rollouts, horizon);
//	}
//	
	//	void sampleSelection() {
	//		MersenneTwisterRNG mersenne_twister;
	//		RandomNumberGenerator* rng = (RandomNumberGenerator*) &mersenne_twister;
	//		rng->manualSeed(323456789);
	//		AbstractPolicy<Vector, int>* policy = new RandomPolicy(environment->getNActions(), rng);
	//		rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> >(environment->getState(), policy, environment, gamma, true);
	//		rollout->UniformSampling(3000, 1);
	////		rollout->Sampling(100, 400);
	//	}
	
	/// Basis function construction based on the cover trees.
	Vector BasisConstruction(const S& state)
	{
		Vector phi(dim);
		if(!strcmp(environment->Name(), "Bike")){
			real psi_hat;
			if(state[4] >= 0) {
				psi_hat = M_PI - state[4];
			}
			else {
				psi_hat = -M_PI - state[4];
			}
			phi.Resize(20);
			phi[0] = 1.0;
			phi[1] = state[2];
			phi[2] = state[3];
			phi[3] = state[2] * state[2];
			phi[4] = state[3] * state[3];
			phi[5] = state[2] * state[3];
			phi[6] = state[0];
			phi[7] = state[1];
			phi[8] = state[0] * state[0];
			phi[9] = state[1] * state[1];
			phi[10] = state[0] * state[1];
			phi[11] = state[2] * state[0];
			phi[12] = state[2] * state[0] * state[0];
			phi[13] = state[2] * state[2] * state[0];
			phi[14] = state[4];
			phi[15] = state[4] * state[4];
			phi[16] = state[4] * state[0];
			phi[17] = psi_hat;
			phi[18] = psi_hat * psi_hat;
			phi[19] = psi_hat * state[0];
		}
		else if(RBFs == NULL && tilings == 0) 
		{
			//		printf("Dimensionnnnnn = %d\n",(int)dim);
			int dim_current = 0;
			for(int i = 0; i<n_actions; ++i) {
				//			printf("Action = %d\n",i);
				//			const CoverTree::Node* node = cover[i]->SelectedNearestNeighbour(state);
				//			int index = dim_current + (node->GetIndex() - 1);
				//			phi(index) = 1.0;
				std::vector< std::pair<int, real> > path = cover[i]->ExternalBasisCreation(state);
				for(uint j = 0; j < path.size(); ++j) {
					//				printf("Index = %d and weight = %f\n",path[j].first,path[j].second);
					int index = dim_current + (path[j].first - 1);
					phi(index) = path[j].second;
				}
				dim_current += cover[i]->GetNumBasisNodes();
//				dim_current += cover[i]->GetNumSamplingNodes();
//				dim_current += cover[i]->GetNumNodes();
			}
		}
		else {
			if(tilings == 0) {
				RBFs->Evaluate(state);
				phi = RBFs->F();
				phi.Resize(dim);
				phi[dim-1] = 1.0; 
			}
			else if(tilings > 0) {
				float scale_div = 1.0 / (float)grids;
//				collision_table ct(dim,1);
				float float_array[2];
				float_array[0] = (float) state[0] / scale_div;
				float_array[1] = (float) state[1] / scale_div;
				
				int tiles_array[tilings];
				GetTiles(tiles_array,tilings,dim,float_array,2);
				for(int i = 0; i<tilings; ++i) {
//					printf("Tiles = %d\n",(int)tiles_array[i]);
					phi[tiles_array[i]] = 1.0;
				}
			}
		}
//		phi.print(stdout);
		return phi;
	}
	
	void Reset() {
//		printf("Dimensions = %d\n",(int)cover[i]->GetNumSamplingNodes());
		if(RBFs==NULL) {
			dim = 0;
			for( int i=0; i<n_actions; ++i) {
//				printf("Dimensions = %d\n",(int)cover[i]->GetNumSamplingNodes());
				dim += cover[i]->GetNumBasisNodes();
//				dim += cover[i]->GetNumSamplingNodes();
//				dim += cover[i]->GetNumNodes();
			}
		}
//		printf("Dimensions = %d\n",dim);
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif

