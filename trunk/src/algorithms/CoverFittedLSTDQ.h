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

#ifndef COVER_FITTED_LSTD_VALUE_ITERATION_Q_H
#define COVER_FITTED_LSTD_VALUE_ITERATION_Q_H

#include "BayesianMultivariateRegression.h"
#include "BasisSet.h"
#include "Environment.h"
#include "Matrix.h"
#include "Random.h"
#include "Rollout.h"
#include "MersenneTwister.h"
#include "RandomPolicy.h"
#include "ContinuousPolicy.h"

///*Fitted Value Iteration Algorithm*/
template <class S, class A>
class CoverFittedLSTDQ {
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
	bool update_samples;		///< Take new random training samples 
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
public:
	CoverFittedLSTDQ(const real& gamma_, const int& N_, const int& M_, Environment<S,A>* environment_, std::vector<CoverTree*> cover_, const int& grids_  = 0, real scale_ = 1.0, bool update_samples_ = false):
	gamma(gamma_),
	N(N_),
	M(M_),
	environment(environment_),
	cover(cover_),
	grids(grids_),
	scale(scale_),
	update_samples(update_samples_)
	{
		n_actions = environment->getNActions();

		if(grids == 0) {
			RBFs = NULL;
			dim = 0;
			for( int i=0; i<n_actions; ++i) {
				dim += cover[i]->GetNumBasisNodes();
				//				dim += cover[i]->GetNumSamplingNodes();
				//				dim += cover[i]->GetNumNodes();
			}
		}
		else {
			//dim	= pow(grids,2.0) + 1;	
			dim = n_actions*(pow(grids,2.0) + 1);
			Vector S_L	= environment->StateLowerBound();
			Vector S_U	= environment->StateUpperBound();
			EvenGrid Discretisation(S_L,S_U,grids);
			RBFs = new RBFBasisSet(Discretisation,scale); 
		}
		//		lambda = 0.01;
		lambda = 0.01;		
		
		weights = Vector(dim);
	}
	const void Update(real threshold = 0.0001, int max_iter = -1) {
		
		Reset();
		
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
		weights = Vector::Unity(dim);
		Vector V(n_actions);
		int n_iter = 0;
		do {
			AA = Matrix::Null(dim,dim);
			b.Clear();
			distance = 0.1;
			real steps = 0.0;
			for(int i=0; i<N; ++i)
			{	
				//Find best action - policy.
				for(a = 0; a < n_actions; a++) {
					environment->Reset();
					environment->setState(states[i]);
				
					phi = BasisConstruction(states[i],a);
					//				phi.print(stdout);
					final = environment->Act(a);
					r = environment->getReward();
				
					if(final) {
						int a_next;
						Vector next_state = cover[a]->GenerateState(states[i]);
						for(a_next = 0; a_next < n_actions; a_next++) { 
							V[a_next] = getValue(next_state,a_next);
						}
						a_next	 = ArgMax(V);
						phi_ = BasisConstruction(next_state, a_next);
						dif	 = phi - (phi_*gamma);
					}
					else {
						dif = phi;
					}	
					steps++;
					Matrix res = OuterProduct(phi,dif);
					AA += res;					
					b  += phi*r;
				}
			}	
			weights = ((1.0/steps)*AA + (lambda*steps)*Matrix::Unity(dim,dim)).Inverse_LU()*(b*(1/steps));
			//	printf("adsfafdsafsdadsfafsdasdf\n");
			weights.print(stdout);
			distance = (pW - weights).L2Norm();
			pW = weights;
			if(max_iter>0) 
				max_iter--;
			n_iter++;
			
			if(update_samples==true) {
				sampleSelection();
			}
		}while(distance > threshold && max_iter != 0);
		printf("#ValueIteration:: ComputeStateValues Exiting at d:%f, n:%d\n", distance, n_iter);
	}
	
	real getValue(const S& state) 
	{
		Vector V(n_actions);
		for(int i = 0; i < n_actions; ++i)
			V[i] = getValue(state, i);
		weights.print(stdout);
		V.print(stdout);
		return Max(V);
	}
	
	real getValue(const S& state, const A& action) 
	{
		Vector phi = BasisConstruction(state,action);
		return Product(weights,phi);
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
	//		MersenneTwisterRNG mersenne_twister;
	//		RandomNumberGenerator* rng = (RandomNumberGenerator*) &mersenne_twister;
	//		rng->manualSeed(323456789);
	//		AbstractPolicy<Vector, int>* policy = new RandomPolicy(environment->getNActions(), rng);
	//		rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> >(environment->getState(), policy, environment, gamma, true);
	//		rollout->UniformSampling(3000, 1);
	////		rollout->Sampling(100, 400);
	//	}
	
	/// Basis function construction based on the cover trees.
	Vector BasisConstruction(const S& state, const A& action)
	{
		assert(action >= 0);
		Vector phi(dim);
		if(RBFs == NULL) 
		{
			int dim_current = 0;
			for(int i = 0; i < action; ++i) {
				dim_current += cover[i]->GetNumBasisNodes();
			}
			
			std::vector< std::pair<int, real> > path = cover[action]->ExternalBasisCreation(state);
			for(uint j = 0; j < path.size(); ++j) {
				int index	= dim_current + (path[j].first - 1);
				phi[index]	= path[j].second; 
			}
		}
		else {
			RBFs->Evaluate(state);
			Vector Phi_state = RBFs->F();
			phi[(RBFs->size() + 1)*action] = 1.0;
			
			for(int i = 0; i<RBFs->size(); ++i)
			{
				phi[(RBFs->size() + 1)*action + i + 1] = Phi_state[i];
			}
			phi[dim-1] = 1.0; 
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

