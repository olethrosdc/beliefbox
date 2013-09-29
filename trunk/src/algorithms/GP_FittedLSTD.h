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

#ifndef GP_FITTED_LSTD_VALUE_ITERATION_H
#define GP_FITTED_LSTD_VALUE_ITERATION_H

#include "GaussianProcess.h"
#include "SparseGaussianProcess.h"
#include "SparseGreedyGaussianProcess.h"
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
class GP_FittedLSTD {
protected:
    real gamma;					///< discount factor
	int N;						///< Number of basis points (sampled states).
	int dim;					///< Dimension of the basis functions
	int n_states;				///< Number of state space dimensions
	int n_actions;				///< Number of actions (Note: what to do for continuous?)
	Vector weights;				///< Model parameters
	Matrix PHI;
	real lambda;				///< Regularization factor
	Matrix pseudo_inv;
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<std::vector<SparseGaussianProcess*> > regression_t;
	std::vector<Matrix> Kernel;
	RBFBasisSet* RBFs;		
	bool update_samples;		///< Take new random training samples 
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout;
public:
	GP_FittedLSTD(const real& gamma_, const int& N_,  Environment<S,A>* environment_, std::vector<std::vector<SparseGaussianProcess*> > regression_t_, RBFBasisSet* RBFs_, bool update_samples_ = false)
		:gamma(gamma_),
		N(N_),
		environment(environment_),
		regression_t(regression_t_),
		RBFs(RBFs_),
		update_samples(update_samples_)
	{
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		
		lambda = 0.01;
	
		n_actions = environment->getNActions();
		n_states  = environment->getNStates();
		
		Kernel = std::vector<Matrix>(n_actions);
			
		dim	= RBFs->size() + 1;	
		
		weights = Vector(dim);
		sampleSelection();
	}
	const void Update(real threshold = 0.00001, int max_iter = 20) {
		
		int a;
	    Matrix AA;
		Vector b(dim);
		Vector dif;
		Vector phi;
		real r;
		bool final;
		real distance;
		Vector T(N);
		Vector pW(dim);
		weights = Vector::Null(dim);
		Vector V(n_actions);
		int n_iter = 0;
//		std::cout << "######### Fitted LSTD #########" << std::endl;
//		std::cout << "####### Iteration = " << n_iter <<  "#######" << std::endl;

//		updateKernels();
		
		do {
//			std::cout << "### Iteration = " << n_iter << std::endl;
			AA = Matrix::Unity(dim,dim)*0;
			b.Clear();
			
			distance = 0.1;
			real steps = 0.0;
			for(int i=0; i<N; ++i)
			{	
				//Find best action - policy.
				for(a = 0; a < n_actions; a++) {
					environment->Reset();
					environment->setState(states[i]);

					final = environment->Act(a);
					r = environment->getReward();
					Vector next_state(n_states);
					for(int d = 0; d<n_states; ++d) {
						next_state(d) = regression_t[a][d]->GeneratePrediction(states[i]);
					}
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
					Vector next_state(n_states);
					for(int d = 0; d<n_states; ++d) {
						next_state(d) = regression_t[a][d]->GeneratePrediction(states[i]);
					}
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
			weights = ((1.0/steps)*AA + (lambda*steps)*Matrix::Unity(dim,dim)).Inverse()*(b*(1/steps));
			distance = (pW - weights).L2Norm();
			pW = weights;
			if(max_iter>0) 
				max_iter--;
			n_iter++;
			
			if(update_samples==true) {
				sampleSelection();
			}
//			printf("#ValueIteration (Threshold = %f) :: ComputeStateValues Exiting at d:%f, n:%d\n", threshold, distance, n_iter);
		}while(distance > threshold && max_iter != 0);
//		printf("#ValueIteration::ComputeStateValues Exiting at d:%f, n:%d\n", distance, n_iter);
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
		
		Vector next_state(n_states);
		for(int d = 0; d<n_states; ++d) {
			next_state(d) = regression_t[action][d]->GeneratePrediction(state);
		}
		temp_v = getValue(next_state);
		
		return (r + gamma*temp_v);
	}
	
	int Act(const S& state)
	{
		Vector Q(n_actions);
		int action = 0;
		real max = getValue(state,0);
		Q(0) = max;

		for(int i = 1; i<n_actions; ++i) {
			Q(i) = getValue(state, i);
			if(Q(i) > max) {
				max = Q(i);
				action = i;
			}
		}
		return action;
	}
	
	void sampleSelection() {
		
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		for(int i=0; i<N; ++i) {
			states.push_back(urandom(S_L, S_U));
		}
	}
	
	void updateKernels() {
		for(int a = 0; a<n_actions; ++a) {
			Kernel[a] = Matrix::Null(N, regression_t[a][0]->getNSamples());
		}
		for(int i=0; i<N; ++i) {
			for(int a = 0; a<n_actions; ++a) {
				Vector k = regression_t[a][0]->Kernel(states[i]);
				Kernel[a].setRow(i, k);
			}
		}
	}
	
	void Reset() {
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif

