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
	FittedValueIteration(const real& gamma_, const int& N_, const int& M_, const int& grids_, Environment<S,A>* environment_, std::vector<BayesianMultivariateRegression*> regression_t_, RBFBasisSet* RBFs_, real scale_ = 1.0, bool update_samples_ = true):
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

		EvenGrid Discretisation(S_L, S_U, grids);
		
		RBFs = new RBFBasisSet(Discretisation,scale); 

		n_actions = environment->getNActions();

		dim			= pow(grids,2.0) + 1;	
		if(RBFs_model!=NULL) {
			dim_model = RBFs_model->size() + 1;
		}
		else {
			dim_model = environment->getNStates() + 1;
		}

		weights = Vector(dim);
		sampleSelection();
	}
	const void Update(real threshold = 0.001, int max_iter = -1) {
		Vector V(n_actions);
		real distance;
		Vector T(N);
		Vector pW(dim);
		weights = pW;
		int n_iter = 0;
		do {
			distance = 0.1;
			
			for(int i=0; i<N; ++i)
			{
				for(int a=0; a<n_actions; ++a) {
					V[a] = 0.0;

					for(int j=0; j<M; ++j) {
						environment->Reset();
						environment->setState(states[i]);
						bool final = environment->Act(a);
						real r = environment->getReward();
						
						if(final) {
							Vector phi;
							if(RBFs_model != NULL) {
								RBFs_model->Evaluate(states[i]);
								phi = RBFs_model->F();
							}
							else {
								phi = states[i];
							}
							phi.Resize(dim_model);
							phi[dim_model-1] = 1.0;
							
							Vector next_state = regression_t[a]->generate(phi);
							real temp_v = getValue(next_state);
							V[a] = V[a] + r + gamma*temp_v;
						}
						else {
							V[a] = V[a] + r;
						}
					}
					V[a] = V[a]/(real)M;
				}
				T[i] = Max(V);
			}
		//	printf("output\n");
//			T.print(stdout);
//			Vector v_old = Transpose(PHI)*weights;
			weights = pseudo_inv*T;
		//	printf("weights\n");
//			weights.print(stdout);
//			printf("PREDICTION\n");
//			Vector VV = Transpose(PHI)*weights;
//			VV.print(stdout);
//			printf("error = %f\n",(abs(T-VV)).Sum());
//			Delta = Max(abs(weights - pW));
			distance = (pW - weights).L2Norm();
			if(max_iter>0) 
				max_iter--;
			n_iter++;
			
			if(update_samples==true) {
				sampleSelection();
			}
		} while(distance > threshold && max_iter != 0);
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
		bool final =environment->Act(action);
		real r = environment->getReward();
		environment->setState(state);

		real temp_v = 0.0;
	
		if(final) {
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
		}
		else {
			temp_v = 0;
		}
		return (r + gamma*temp_v);
	}
	
	void sampleSelection() {
		real lambda = 0.001; //regularization factor

		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		PHI = Matrix(dim,N);
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
			RBFs->Evaluate(state);
			Vector phi = RBFs->F();
			phi.Resize(dim);
			phi[dim-1] = 1.0;
			PHI.setColumn(i,phi);
		}
		pseudo_inv = (PHI*Transpose(PHI) + lambda*Matrix::Unity(dim,dim)).Inverse()*PHI;
	}
	
	void Reset() {
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif
