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

#ifndef FITTED_LSTD_VALUE_ITERATION_H
#define FITTED_LSTD_VALUE_ITERATION_H

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
class FittedLSTD {
protected:
    real gamma;					///< discount factor
	int N;						///< Number of basis points
	int M;						///< Number of sampled next states
	int grids;					///< Number of grids of basis function
	int dim;					///< Dimension of the basis functions
	int dim_model;				///< Dimension of the model basis functions
	int n_actions;				///< Number of actions (Note: what to do for continuous?)
	Vector weights;				///< Model parameters
	Matrix PHI;
	real lambda;				///< Regularization factor
	Matrix pseudo_inv;
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<BayesianMultivariateRegression*> regression_t;
	RBFBasisSet* RBFs_model;	///< The Radial basis functions that used in the model.
	RBFBasisSet* RBFs;		
	real scale; ///< scale parameter
	bool update_samples;		///< Take new random training samples 
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
public:
	FittedLSTD(const real& gamma_, const int& N_, const int& M_, const int& grids_, Environment<S,A>* environment_, std::vector<BayesianMultivariateRegression*> regression_t_, RBFBasisSet* RBFs_, real scale_ = 1.0, bool update_samples_ = false):
	gamma(gamma_),
	N(N_),
	M(M_),
	grids(grids_),
	environment(environment_),
	regression_t(regression_t_),
	RBFs_model(RBFs_),
	scale(scale_),
	update_samples(update_samples_)
	{
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		
		lambda = 0.01;
		EvenGrid Discretisation(S_L,S_U,grids);
		RBFs = new RBFBasisSet(Discretisation,scale); 
		n_actions = environment->getNActions();
		
		dim	= pow(grids,environment->getNStates()) + 1;	
		if(RBFs_model!=NULL) {
			dim_model = RBFs_model->size() + 1;
		}
		else {
			dim_model = environment->getNStates() + 1;
		}
		
		weights = Vector(dim);
		sampleSelection();
	}
	const void Update(real threshold = 0.0001, int max_iter = -1) {
		
		int a;
	    Matrix AA;
		Vector b(dim);
		Vector dif;
		Vector phi_model;
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
				Vector phi_model;
				if(RBFs_model != NULL) {
					RBFs_model->Evaluate(states[i]);
					phi_model = RBFs_model->F();
				}
				else {
					phi_model = states[i];
				}
				phi_model.Resize(dim_model);
				phi_model[dim_model-1] = 1.0;
				//Find best action - policy.
				for(a = 0; a < n_actions; a++) {
					environment->Reset();
					environment->setState(states[i]);
					
					final = environment->Act(a);
					r = environment->getReward();
					Vector next_state = regression_t[a]->generate(phi_model);
					V[a] = getValue(next_state) + r;
				}
				
				a = ArgMax(V);
				
				environment->Reset();
				environment->setState(states[i]);
								
				RBFs->Evaluate(states[i]);
				phi = RBFs->F();
				phi.Resize(dim);
				phi[dim-1] = 1.0; 
				
				final = environment->Act(a);
				r = environment->getReward();
				
				if(final) {
					Vector next_state = regression_t[a]->generate(phi_model);
					RBFs->Evaluate(next_state);
					Vector phi_ = RBFs->F();
					phi_.Resize(dim);
					phi_[dim-1] = 1.0;
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
			weights = ((1.0/steps)*AA + (lambda*steps)*Matrix::Unity(dim,dim)).Inverse_LU()*(b*(1/steps));
			
			distance = (pW - weights).L2Norm();
			pW = weights;
			if(max_iter>0) 
				max_iter--;
			n_iter++;
			
			if(update_samples==true) {
				sampleSelection();
			}
		}while(distance > threshold && max_iter != 0);
		printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", distance, n_iter);
	}
	
	real getValue(const S& state) 
	{
		RBFs->Evaluate(state);
		Vector phi = RBFs->F();
		phi.Resize(dim);
		phi[dim-1] = 1.0;
		
		return Product(weights,phi);
	}
	
	real getValue(const S& state, const A& action) 
	{
		environment->Reset();
		environment->setState(state);
		environment->Act(action);
		real r		= environment->getReward();
//		environment->Reset();

		environment->setState(state);
		
		real temp_v = 0.0;
		
		Vector phi;
		if(RBFs_model != NULL) {
			RBFs_model->Evaluate(state);
			phi = RBFs_model->F();
		}
		else {
			phi = state;
		}
		phi.Resize(dim_model);
		phi[dim_model-1] = 1.0;
		Vector next_state = regression_t[action]->generate(phi);
		temp_v = getValue(next_state);
	
		return (r + gamma*temp_v);
	}
	
	void sampleSelection() {
		
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		PHI = Matrix(dim,N);
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
		}
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
	
	void Reset() {
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif

