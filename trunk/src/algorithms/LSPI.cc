/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "LSPI.h"

LSPI::LSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, RBFBasisSet* bfs_, Rollout<Vector,int,AbstractPolicy<Vector, int> >* Samples_)
	:gamma(gamma_), Delta(Delta_), n_dimension(n_dimension_), n_actions(n_actions_), max_iteration(max_iteration_),bfs(bfs_), Samples(Samples_)
{
	assert(gamma>=0 && gamma <=1);
	n_basis = n_actions*bfs->size();
	A.Resize(n_basis, n_basis);
	b.Resize(n_basis);
	w.Resize(n_basis);
	policy = new FixedContinuousPolicy( n_dimension_, n_actions_, bfs_);
}

Vector LSPI::BasisFunction(Vector state, int action)
{
	bfs->Evaluate(state);
	Vector Phi_state = bfs->F();
	Vector Phi(n_basis);
	for(int i = 0; i<bfs->size(); ++i)
	{
		Phi[bfs->size()*action + i] = Phi_state[i];
	}
	return Phi;
}

void LSPI::LSTDQ()
{
	Vector Phi_;
	Vector Phi;
	Matrix res;
    A.Clear();
	b.Clear();

	for(int i=0; i<Samples->getNRollouts(); ++i)
	{
		for(int j=0; j<Samples->getNSamples(i); ++j)
		{
			Phi_ = BasisFunction(Samples->getState(i,j), Samples->getAction(i,j));
			if(Samples->getEndsim(i,j)){
				res = OuterProduct(Phi_, Phi_);
			}
			else{
				Phi = BasisFunction(Samples->getNextState(i,j),policy->SelectAction(Samples->getNextState(i,j)));
				res = OuterProduct(Phi_,(Phi_ - (Phi*gamma)));
			}
			A = A + res;
			b = b + Phi_*Samples->getReward(i,j);
		}
	}
	const Matrix w_ = A.Inverse_LU();
	w = w_*b;
}

void LSPI::PolicyIteration()
{
	Vector old_w;
	real distance;
	int iteration = 0;
	while(1)
	{
//		printf("Policy Evaluation\n");
		LSTDQ();
//		printf("Policy Improvement\n");
		old_w = policy->getWeights();
		policy->Update(w);
		
		iteration++;
		//Stop criterion
		distance = (old_w - w).L2Norm();
		if( distance < Delta || iteration > max_iteration){
			if(distance > Delta){
				printf("LSPI finished after %d iterations without convergence into a fixed point\n",iteration);
			}
			else{
				printf("LSPI converged after %d iterations",iteration);
			}	
			break;
		}
	}
}

real LSPI::getValue(Vector state, int action)
{
	return Product(BasisFunction(state,action),w);
}

