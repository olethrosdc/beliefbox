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
	real lambda;				///< Regularization factor
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<BayesianMultivariateRegression*> regression_t;
	RBFBasisSet* RBFs_model;	///< The Radial basis functions that used in the model.
	RBFBasisSet* RBFs;		
	real scale;					///< scale parameter
	bool update_samples;		///< Take new random training samples 
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
		///Upper and Lower environment's bounds.
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		
		lambda = 0.0001; ///Regularize parameter.
		///Basis function construction for the API model
		EvenGrid Discretisation(S_L,S_U,grids);
		RBFs = new RBFBasisSet(Discretisation,scale); 
		dim	= RBFs->size() + 1;	

		n_actions = environment->getNActions();
		
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
		Vector state, next_state;
		Vector phi_model, phi, phi_;
		real r;
		real distance = 0.1;
		bool final;
		Vector pW(dim); //previous weigths
		weights = Vector::Null(dim);
		Vector V(n_actions);
		int n_iter = 0;
		do {
			AA = Matrix::Unity(dim,dim)*0.0;
			b.Clear();

			for(int i=0; i<N; ++i)
			{	
				phi_model = BasisModelCreation(states[i]);
				
				//Find best action - policy.
				for(a = 0; a < n_actions; a++) {
					environment->Reset();
					environment->setState(states[i]);
					
					final = environment->Act(a);
					r = environment->getReward();
					///We find the next state according to the learned environment model
					if(final) {
						next_state = environment->getState();
						V[a] = getValue(next_state) + r;
					} else {
						V[a] = r;
					}
				}
				
				a = ArgMax(V);
				
				environment->Reset();
				environment->setState(states[i]);
				///In this point we calculate the basis function for the collected state i				
				phi = BasisAPICreation(states[i]);
				
				final = environment->Act(a);
				r = environment->getReward();
				///Least Square Temporal Difference Update
				if(final) {
					///We find the next state according to the learned environment model
					next_state = regression_t[a]->generate(phi_model);
					///In this point we calculate the basis function for the next state			
					phi_ = BasisAPICreation(next_state);
					dif = phi - (phi_*gamma);
				}
				else {
					dif = phi;
				}	
				Matrix res = OuterProduct(phi,dif);
				AA = AA + res;					
				b = b + phi*r;
			}	
			
			weights = ((1.0/N)*AA + (lambda)*Matrix::Unity(dim,dim)).Inverse()*(b*(1.0/N));
			
			distance = (pW - weights).L2Norm();
			pW = weights;
		
			max_iter--;
			n_iter++;
			/// If "true" we collect new samples on each iteration
			if(update_samples==true) {
				sampleSelection();
			}
		}while(distance > threshold && max_iter != 0);
		printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", distance, n_iter);
	}
	
	real getValue(const S& state) 
	{
		Vector phi = BasisAPICreation(state);
		return Product(weights,phi);
	}
	
	real getValue(const S& state, const A& action) 
	{
		bool endsim = environment->getEndsim();
		Vector true_state = environment->getState();
		environment->Reset();
		environment->setState(state);
		bool final = environment->Act(action);
		real r		= environment->getReward();

		environment->Reset();
		environment->setEndsim(endsim);
		environment->setState(true_state);
		
		real temp_v = 0.0;
		if(final) {
			Vector phi = BasisModelCreation(state);
			Vector next_state = regression_t[action]->generate(phi);
			temp_v = getValue(next_state);
		}
		return (r + gamma*temp_v);
	}
	///SampleSelection collects a number of samples, uniformly random
	void sampleSelection() {
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
		}
	}
	// BasisModelCreation returns the basis function for state s that used for the Model prediction algorithm
	Vector BasisModelCreation(const Vector& s) {
		Vector phi;
		if(RBFs_model != NULL) {
			RBFs_model->Evaluate(s);
			phi = RBFs_model->F();
		} else {
			phi = s;
		}
		phi.Resize(dim_model);
		phi[dim_model-1] = 1.0;
		return phi;
	}
	// BasisAPICreation returns the basis function for state s that used for the API algorithm
	Vector BasisAPICreation(const Vector& s) {
		RBFs->Evaluate(s);
		Vector phi = RBFs->F();
		phi.Resize(dim);
		phi[dim-1] = 1.0;
		return phi;
	}
	
	void Reset() {
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif

