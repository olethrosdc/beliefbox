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

#include "OnlineLSPI.h"

OnlineLSPI::OnlineLSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, RBFBasisSet* bfs_)
:gamma(gamma_),
Delta(Delta_), 
n_dimension(n_dimension_), 
n_actions(n_actions_), 
max_iteration(max_iteration_),
bfs(bfs_), 
policy(n_dimension, n_actions, bfs)
{
	assert(gamma>=0 && gamma <=1);
	n_basis = n_actions*(bfs->size() + 1);
	algorithm = 1;
	A.Resize(n_basis, n_basis);
	A = Matrix::Unity(n_basis,n_basis) * 1e-6;
	b.Resize(n_basis);
	w.Resize(n_basis);
}

OnlineLSPI::OnlineLSPI(real gamma_, real Delta_, int n_dimension_, int n_actions_, int max_iteration_, int algorithm_, RBFBasisSet* bfs_)
:gamma(gamma_),
Delta(Delta_), 
n_dimension(n_dimension_), 
n_actions(n_actions_), 
max_iteration(max_iteration_),
algorithm(algorithm_),
bfs(bfs_), 
policy(n_dimension, n_actions, bfs)
{
	assert(gamma>=0 && gamma <=1);
	assert(algorithm>=1 && algorithm<=2);
	n_basis = n_actions*(bfs->size() + 1);
	A.Resize(n_basis, n_basis);
	A = Matrix::Unity(n_basis,n_basis) * 1e-6;
	b.Resize(n_basis);
	w.Resize(n_basis);
}

OnlineLSPI::~OnlineLSPI()
{
}

Vector OnlineLSPI::BasisFunction(const Vector& state, int action)
{
	bfs->Evaluate(state);
	Vector Phi_state = bfs->F();
	Vector Phi(n_basis);
	Phi[(bfs->size() + 1)*action] = 1.0;
	
	for(int i = 0; i<bfs->size(); ++i)
	{
		Phi[(bfs->size() + 1)*action + i + 1] = Phi_state[i];
	}
	return Phi;
}
void OnlineLSPI::LSTD(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update) 
{ 
	if( algorithm == 1 )
		LSTDQ(state, action, reward, state_, action_, endsim, update);
	else if( algorithm == 2 )
		LSTDQ_OPT(state, action, reward, state_, action_, endsim, update);
}
void OnlineLSPI::LSTDQ(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update) 
{
	Vector Phi_;
	Vector Phi;
	Matrix res;
	
	Phi_ = BasisFunction(state, action);
	if(endsim) {
		res = OuterProduct(Phi_, Phi_);
	}
	else {
		Phi = BasisFunction(state_,action_);
		res = OuterProduct(Phi_,(Phi_ - (Phi*gamma)));
	}
	
	A += res;
	b += Phi_*reward;
}
void OnlineLSPI::LSTDQ_OPT(const Vector& state, const int& action, const real& reward, const Vector& state_, const int& action_, const bool& endsim, const bool& update)
{
	Vector Phi_;
	Vector Phi;
	Vector Phi_dif;
	Matrix res;
	real d = 0.000001;
	A = Matrix::Unity(n_basis,n_basis) * (1/d);
	b.Clear();
	
	Phi_ = BasisFunction(state, action);
	if(endsim){
		Phi_dif = Phi_;
	}
	else{
		Phi = BasisFunction(state_, action_);
		Phi_dif = Phi_ - (Phi*gamma);
	}
	res = OuterProduct(Phi_,Phi_dif);
	const Matrix p = A;
	real v = Product(p*Phi_,Phi_dif);
	A -= (((A*res)*A) / (v + 1));
	b += Phi_ * reward;
}
void OnlineLSPI::Update()
{
	if(algorithm == 1) {
		const Matrix w_ = A.Inverse_LU();
		w = w_*b;
	}
	else if( algorithm == 2) {
		const Matrix w_ = A;
		w = w_*b;
	}
	policy.Update(w);
}
void OnlineLSPI::Reset()
{
	A = Matrix::Unity(n_basis,n_basis) * 1e-6;
	b = Vector::Null(n_basis);
	w = Vector::Null(n_basis);	
}
real OnlineLSPI::getValue(const Vector& state, int action)
{
	return Product(BasisFunction(state,action),w);
}

