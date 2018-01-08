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

#ifndef FITTED_VALUE_ITERATION_H
#define FITTED_VALUE_ITERATION_H

#include "BayesianMultivariateRegression.h"
#include "BasisSet.h"
#include "Environment.h"
#include "Matrix.h"

///*Fitted Value Iteration Algorithm*/
template <class S, class A>
class FittedValueIteration {
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
	Matrix PHI;
	Matrix pseudo_inv;
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<BayesianMultivariateRegression*> regression_t;
	RBFBasisSet* RBFs_model;	///< The Radial basis functions that used in the model.
	RBFBasisSet* RBFs;
	real scale;
	bool update_samples;		///< Take new random training samples 
public:
	FittedValueIteration(const real& gamma_,
						 const int& N_,
						 const int& M_,
						 const int& grids_,
						 Environment<S,A>* environment_,
						 std::vector<BayesianMultivariateRegression*> regression_t_,
						 RBFBasisSet* RBFs_,
						 real scale_ = 1.0,
						 bool update_samples_ = false):
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
		
		lambda = 0.01; ///Regularize parameter.
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
	const void Update(real threshold = 0.000001, int max_iter = -1) {
		real distance;
		Vector phi, next_state;
		Vector T(N);
		Vector pW(dim);
		weights = pW;
		int n_iter = 0;
		do {
			distance = 0.1;
			
			for(int i=0; i<N; ++i)
			{
				Vector V(n_actions);
				for(int a=0; a<n_actions; ++a) {
					V[a] = 0.0;

					for(int j=0; j<M; ++j) {
						environment->Reset();
						environment->setState(states[i]);
						bool final = environment->Act(a);
						real r = environment->getReward();
						if(final) {
							phi = BasisModelCreation(states[i]);
							next_state = regression_t[a]->generate(phi);

							real temp_v = getValue(next_state);
							V[a] = V[a] + (r + gamma*temp_v);
						}
						else {
							V[a] = V[a] + r;
						}
					}
					V[a] = V[a]/(real)M;
				}
				T[i] = Max(V);
			}
			//Update the weights parameters.
			weights = pseudo_inv*T;
			//Find the distance between the weights of successive learning iterations.
			distance = (pW - weights).L2Norm();
			pW = weights;
			max_iter--;
			n_iter++;
			/// If "true" we collect new samples on each iteration
			if(update_samples==true) {
				sampleSelection();
			}
		} while(distance > threshold && max_iter != 0);
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
		bool final =environment->Act(action);
		real r = environment->getReward();
		
		environment->Reset();
		environment->setEndsim(endsim);
		environment->setState(true_state);
		
		real temp_v = 0.0;
	
		if(final) {
			Vector phi			= BasisModelCreation(state);
			Vector next_state	= regression_t[action]->generate(phi);
			temp_v = getValue(next_state);
		}
		else {
			temp_v = 0;
		}
		return (r + gamma*temp_v);
	}
	///SampleSelection collects a number of samples, uniformly random
	void sampleSelection() {
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		PHI = Matrix(dim,N);
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
			Vector phi = BasisAPICreation(state);
			PHI.setColumn(i,phi);
		}
		pseudo_inv = (PHI*Transpose(PHI) + lambda*Matrix::Unity(dim,dim)).Inverse()*PHI;
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
