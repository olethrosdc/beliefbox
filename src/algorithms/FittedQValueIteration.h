/* -*- Mode: C++; -*- */
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef FITTED_Q_VALUE_ITERATION_H
#define FITTED_Q_VALUE_ITERATION_H

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
class FittedQValueIteration {
protected:
	real gamma;					///< discount factor
	int N;						///< Number of basis points
	int M;						///< Number of sampled next states
	int grids;					///< Number of grids of basis function
	int dim;					///< Dimension of the basis functions
	int dim_model;				///< Dimension of the model basis functions
	int n_actions;				///< Number of actions (Note: what to do for continuous?)
	Vector weights;///< Model parameters
	Matrix PHI;
	Matrix pseudo_inv;
	Environment<S, A>* environment;
	std::vector<BayesianMultivariateRegression*> regression_t;
	RBFBasisSet* RBFs_model;	///< The Radial basis functions that used in the model.
	RBFBasisSet* RBFs;		
	real scale;
	bool update_samples;		///< Take new random training samples 
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
	std::vector<Vector> states;	///< set of representative states
public:
	FittedQValueIteration(const real& gamma_,
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
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		
		EvenGrid Discretisation(S_L,S_U,grids);
		RBFs		= new RBFBasisSet(Discretisation,scale); 
		n_actions	= environment->getNActions();
		
		dim	= n_actions*(RBFs->size() + 1);	
		
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

		Vector T(rollout->getNSamples());
		Vector pW(dim);
		printf("Number of states = %d\n",rollout->getNSamples());
		weights = pW;
		int a;
		int counter;
		Vector state;
		int n_iter = 0;
		do {
			counter = 0;
			distance = 0.1;
			for(int i=0; i<rollout->getNRollouts(); ++i)
			{
				for(int j=0; j<rollout->getNSamples(i); ++j) {
					a		= rollout->getAction(i,j);
					state	= rollout->getState(i,j);
				    V[a]	= 0.0;
					for(int j=0; j<M; ++j) {
						environment->Reset();
						environment->setState(state);
						bool final	= environment->Act(a);
						real r		= environment->getReward();
						
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
							
							Vector next_state = regression_t[a]->generate(phi);
							real temp_v = getValue(next_state);
//							printf("Temp => %f\n",temp_v);
							V[a] = V[a] + r + gamma*temp_v;
						}
						else {
							V[a] = V[a] + r;
						}
					}
					counter++;
					V[a] = V[a]/(real)M;
					T[counter] = Max(V);
				}
			}
//			printf("output\n");
//			T.print(stdout);
//			weights.print(stdout);
			pW = weights;
//			Vector v_old = Transpose(PHI)*weights;
//			printf("weights\n");
			weights = pseudo_inv*T;
//			weights.print(stdout);
//			printf("PREDICTION\n");
			Vector VV = Transpose(PHI)*weights;
			
//			VV.print(stdout);
//			printf("error = %f\n",(abs(T-VV)).Sum());
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
	Vector BasisFunction(const S& state, const A& action)
	{
		RBFs->Evaluate(state);
		Vector Phi_state = RBFs->F();
		Vector Phi(dim);
		Phi[(RBFs->size() + 1)*action] = 1.0;
		
		for(int i = 0; i<RBFs->size(); ++i)
        {
            Phi[(RBFs->size() + 1)*action + i + 1] = Phi_state[i];
        }
		return Phi;
		
	//	RBFs->Evaluate(state);
//		Vector Phi = RBFs->F();
//
//		Phi.Resize(dim_model);
//		Phi[dim-1] = 1.0;
//		
//		return Phi;
	}
	real getValue(const S& state)
	{
		Vector Q(n_actions);
		for(int i = 0; i<n_actions; i++) 
			Q[n_actions] = getValue(state,i);
		return Max(Q);
	}
	real getValue(const S& state, const A& action)
	{
		Vector phi = BasisFunction(state,action);
		return Product(weights,phi);
	}
	//void sampleSelection() {
//		real lambda = 1; //regularization factor
//			
//		Vector S_L	= environment->StateLowerBound();
//		Vector S_U	= environment->StateUpperBound();
//		states.clear();
//		for(int i=0; i<N; ++i) 
//			states.push_back(urandom(S_L, S_U));
//	}
	void sampleSelection() {
		real lambda =  1000; //regularization factor

		MersenneTwisterRNG mersenne_twister;
		RandomNumberGenerator* rng = (RandomNumberGenerator*) &mersenne_twister;
		rng->manualSeed(323456789);
		AbstractPolicy<Vector, int>* policy = new RandomPolicy(environment->getNActions(), rng);
		rollout = new Rollout<Vector, int, AbstractPolicy<Vector, int> >(environment->getState(), policy, environment, gamma, true);
		rollout->UniformSampling(N, 1);
		//		rollout->Sampling(100, 400);
		int count = 0;
		PHI = Matrix(dim,rollout->getNSamples());
		lambda *= rollout->getNSamples();
		//		Vector counter(3);
		for(int i=0; i<rollout->getNRollouts(); ++i)
		{
			for(int j=0; j<rollout->getNSamples(i); ++j) {
				Vector phi = BasisFunction(rollout->getState(i,j),rollout->getAction(i,j));
				PHI.setColumn(count,phi);
				count++;
			}
		}
//		printf("Matrix\n");
		
//		(PHI*Transpose(PHI)).print(stdout);
		pseudo_inv = (PHI*Transpose(PHI) + lambda*Matrix::Unity(dim,dim)).Inverse_LU()*PHI;
	}
	void Reset() {
		sampleSelection();
		weights = Vector(dim);
	}		
};
#endif

