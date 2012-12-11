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

#include "LSTDQ.h"

LSTDQ::LSTDQ(real gamma_,
             int n_dimension_,
             int n_actions_,
             RBFBasisSet& bfs_,
             Demonstrations<Vector, int>& Samples_)
	:gamma(gamma_),
	 n_dimension(n_dimension_), 
	 n_actions(n_actions_), 
	 bfs(bfs_), Samples(Samples_), 
	 policy(n_dimension, n_actions, &bfs)
{
    assert(gamma>=0 && gamma <=1);
    n_basis = n_actions*(bfs.size() + 1);
    algorithm = 1;
    A.Resize(n_basis, n_basis);
    b.Resize(n_basis);
    w.Resize(n_basis);
}

LSTDQ::LSTDQ(real gamma_,
             int n_dimension_,
             int n_actions_,
             int algorithm_,
             RBFBasisSet& bfs_,
             Demonstrations<Vector, int>& Samples_)
	:gamma(gamma_),
	 n_dimension(n_dimension_), 
	 n_actions(n_actions_), 
	 algorithm(algorithm_),
	 bfs(bfs_), Samples(Samples_), 
	 policy(n_dimension, n_actions, &bfs)
{
    assert(gamma>=0 && gamma <=1);
    assert(algorithm>=1 && algorithm<=2);
    n_basis = n_actions*(bfs.size() + 1);
    A.Resize(n_basis, n_basis);
    b.Resize(n_basis);
    w.Resize(n_basis);
}

LSTDQ::~LSTDQ()
{
}

Vector LSTDQ::BasisFunction(const Vector& state, int action) const
{
    bfs.Evaluate(state);
    Vector Phi_state = bfs.F();
    Vector Phi(n_basis);
    Phi[(bfs.size() + 1)*action] = 1.0;
    
    for(int i = 0; i<bfs.size(); ++i)
        {
            Phi[(bfs.size() + 1)*action + i + 1] = Phi_state[i];
        }
    return Phi;
}

void LSTDQ::Calculate()
{
    Vector Phi_;
    Vector Phi;
    Matrix res;
    A.Clear();
    b.Clear();

    for(uint i=0; i<Samples.size(); ++i) {
        for(uint t=0; t<Samples.length(i) - 1; ++t) {
            Phi_ = BasisFunction(Samples.state(i,t), Samples.action(i,t));
            Vector s = Samples.state(i, t+1);
            Phi = BasisFunction(s, policy.SelectAction(s));
            res = OuterProduct(Phi_,(Phi_ - (Phi*gamma)));
            A += res;
            b += Phi_*Samples.reward(i,t);
        }
    }
    const Matrix w_ = A.Inverse_LU();
    w = w_*b;
}
void LSTDQ::Calculate_Opt()
{
    Vector Phi_;
    Vector Phi;
    Vector Phi_dif;
    Matrix res;
    real d = 0.000001;
    A = Matrix::Unity(n_basis,n_basis) * (1/d);
    b.Clear();
    
    for(uint i=0; i<Samples.size(); ++i) {
        for(uint t=0; t<Samples.length(i); ++t) {
            Phi_ = BasisFunction(Samples.state(i,t), Samples.action(i,t));
            Vector s = Samples.state(i, t+1);
            Phi = BasisFunction(s,policy.SelectAction(s));
            Phi_dif = Phi_ - (Phi*gamma);
            res = OuterProduct(Phi_,Phi_dif);
            const Matrix p = A;
            real v = Product(p*Phi_,Phi_dif);
            A -= (((A*res)*A) / (v + 1));
            b += Phi_ * Samples.reward(i,t);
        }
    }
    const Matrix w_ = A;
    w = w_*b;
}

real LSTDQ::getValue(const Vector& state, int action) const
{
    return Product(BasisFunction(state,action),w);
}

