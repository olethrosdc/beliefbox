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

#ifndef FITTED_VALUE_ITERATION_H
#define FITTED_VALUE_ITERATION_H

#include "CoverTree.h"
#include "BasisSet.h"
#include "Environment.h"
#include "Matrix.h"

///*Fitted Value Iteration Algorithm*/
template <class S, class A>
class CoverFittedValueIteration {
protected:
    real gamma;					///< discount factor
	int N;						///< Number of basis points
	int M;						///< Number of sampled next states
	int grids;					///< Number of grids of basis function
	int dim;					///< Dimension of the basis functions
	int n_actions;				///< Number of actions (Note: what to do for continuous?)
	Vector weights;				///< Model parameters
	Matrix PHI;
	Matrix pseudo_inv;
    std::vector<Vector> states;	///< set of representative states
	Environment<S, A>* environment;
	std::vector<CoverTree*> cover;
	bool update_samples;		///< Take new random training samples 
public:
	CoverFittedValueIteration(const real& gamma_, const int& N_, const int& M_, Environment<S,A>* environment_, std::vector<CoverTree*> cover_, bool update_samples_ = false):
	gamma(gamma_),
	N(N_),
	M(M_),
	environment(environment_),
	cover(cover_),
	update_samples(update_samples_)
	{		
		n_actions = environment->getNActions();
		printf("adfaf\n");
		dim = 0;
		for( int i=0; i<n_actions; ++i) {
			dim += cover[i]->GetNumSamplingNodes();
		}
		weights = Vector(dim);
		printf("adfaf dimension = %d\n",dim);
		printf("YEAH\n");

	}
	const void Update(real threshold = 0.001, int max_iter = -1) {
//		sampleSelection();
		Reset();
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
							Vector next_state	= cover[a]->GenerateState(states[i]);
							real temp_v			= getValue(next_state);
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
			printf("output\n");
			T.print(stdout);
			Vector v_old = Transpose(PHI)*weights;
			weights = pseudo_inv*T;
			printf("weights\n");
			weights.print(stdout);
			printf("PREDICTION\n");
			Vector VV = Transpose(PHI)*weights;
			VV.print(stdout);
			printf("error = %f\n",(abs(T-VV)).Sum());
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
		Vector phi = BasisConstruction(state);
		return Product(weights,phi);
	}
	
	real getValue(const S& state, const A& action) 
	{
		environment->Reset();
		environment->setState(state);
		bool final	=environment->Act(action);
		real r		= environment->getReward();
		environment->setState(state);
		
		real temp_v = 0.0;
		
		if(final) {
			Vector next_state = cover[action]->GenerateState(state);
			temp_v = getValue(next_state);
		}
		else {
			temp_v = 0;
		}
		return (r + gamma*temp_v);
	}
	void sampleSelection() {
		real lambda = 0.01; //regularization factor
		
		Vector S_L	= environment->StateLowerBound();
		Vector S_U	= environment->StateUpperBound();
		states.clear();
		PHI = Matrix(dim,N);
		for(int i=0; i<N; ++i) {
			Vector state = urandom(S_L, S_U);
			states.push_back(state);
			Vector phi = BasisConstruction(state);			
			PHI.setColumn(i,phi);
		}
		pseudo_inv = (PHI*Transpose(PHI) + lambda*Matrix::Unity(dim,dim)).Inverse()*PHI;
	}
	/// Basis function construction based on the cover trees.
	Vector BasisConstruction(const S& state)
	{
		Vector phi(dim);
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
			dim_current += cover[i]->GetNumSamplingNodes();
		}
		return phi;
	}
	
	void Reset() {
		dim = 0;
		for( int i=0; i<n_actions; ++i) {
//			printf("Dimensions = %d\n",(int)cover[i]->GetNumSamplingNodes());
			dim += cover[i]->GetNumSamplingNodes();
		} 
//		printf("Dimensions = %d\n",dim);
		sampleSelection();
		weights = Vector(dim);
	}
};

#endif
