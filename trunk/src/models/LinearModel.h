/// -*- Mode: c++ -*-
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef LINEAR_MODEL_H
#define LINEAR_MODEL_H

#include "BayesianMultivariateRegression.h"
#include "BasisSet.h"
#include "Environment.h"
#include "MultivariateNormal.h"
#include "Matrix.h"

template <typename S, typename A>
class LinearModel 
{		
protected:
	int m;					///< Input dimensions
	int d;					///< Output dimensions of the transition model
	std::vector<Matrix> M;	///< Mean Matrix of the transition model
	std::vector<Matrix> V;	///< Covariance Matrix of the transition model
	std::vector<Matrix> SP;  
	std::vector<Vector> MR; ///< Mean of the reward model
	std::vector<Matrix> VR; ///< Covariance Matrix of the reward model
	std::vector<Matrix> SR;
	RBFBasisSet* RBFs;
	Environment<S, A>* environment;
	std::vector<BayesianMultivariateRegression*> regression_t;
public:
	LinearModel(const int& m_, const int& d_, RBFBasisSet* RBFs_, Environment<S,A>* environment_, std::vector<BayesianMultivariateRegression*> regression_t_) : RBFs(RBFs_), environment(environment_), regression_t(regression_t_)
	{
		m = m_;
		d = d_;
		int n_actions = environment->getNActions();
		M.resize(n_actions);
		V.resize(n_actions);
		SP.resize(n_actions);
		MR.resize(n_actions);
		VR.resize(n_actions);
		SR.resize(n_actions);
		for(int i=0; i<n_actions; ++i) {
			M[i]	= Matrix::Null(d,m);
			V[i]	= Matrix::Unity(d,d);
			SP[i]	= Matrix::Unity(m,m);
			MR[i]	= Vector::Null(m);
			VR[i]	= Matrix::Unity(1,1);
			SR[i]	= Matrix::Unity(1,1);
		}
	}
	LinearModel(const int& m_, const int& d_, Environment<S,A>* environment_) : environment(environment_)
	{
		m = m_;
		d = d_;
		int n_actions = environment->getNActions();
		M.resize(n_actions);
		V.resize(n_actions);
		SP.resize(n_actions);
		MR.resize(n_actions);
		VR.resize(n_actions);
		SR.resize(n_actions);
		for(int i=0; i<n_actions; ++i) {
			M[i]	= Matrix::Null(d,m);
			V[i]	= Matrix::Unity(d,d);
			SP[i]	= Matrix::Unity(m,m);
			MR[i]	= Vector::Null(m);
			VR[i]	= Matrix::Unity(1,1);
			SR[i]	= Matrix::Unity(1,1);
		}
		RBFs = NULL;
	}
	LinearModel(const std::vector<Matrix>& M_, 
				const std::vector<Matrix>& V_, 
				const std::vector<Matrix>& SP_,
				const std::vector<Vector>& MR_, 
				const std::vector<Matrix>& VR_,
				const std::vector<Matrix>& SR_,
				RBFBasisSet* RBFs_, 
				Environment<S,A>* environment_)
		: M(M_), V(V_), SP(SP_), MR(MR_), VR(VR_), SR(SR_), RBFs(RBFs_), environment(environment_)
	{
		d = M[0].Rows();
		m = M[0].Columns();
		//*TODO*/
	}
	LinearModel(const std::vector<Matrix>& M_, 
				const std::vector<Matrix>& V_, 
				const std::vector<Matrix>& SP_,
				const std::vector<Vector>& MR_, 
				const std::vector<Matrix>& VR_,
				const std::vector<Matrix>& SR_,
				Environment<S,A>* environment_)
	: M(M_), V(V_), SP(SP_), MR(MR_), VR(VR_), SR(SR_), environment(environment_)
	{
		d = M[0].Rows();
		m = M[0].Columns();
		RBFs = NULL;
		//*TODO*/
	}
	void Reset()
	{
		int n_actions = environment->getNActions();
		for(int i=0; i<n_actions; ++i) {
			M[i]	= Matrix::Null(d,m);
			V[i]	= Matrix::Unity(d,d);
			SP[i]	= Matrix::Unity(m,m);
			MR[i]	= Vector::Null(m);
			VR[i]	= Matrix::Unity(1,1);
			SR[i]	= Matrix::Unity(1,1);
		}
	}
	const void SetStatePredictionMean(const Matrix& M_, const int& a)
	{
		assert(M_.Rows() == d && M_.Columns() == m);
		M[a] = M_;
	}
	const void SetStatePredictionMean(const std::vector<Matrix>& M_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert( M_[i].Rows() == d && M_[i].Columns() == m);
		}
		M = M_;
	}
	Matrix GetStatePredictionVar(const int& a)
	{
		return V[a];
	}
	const void SetStatePredictionVar(const Matrix& V_, const int& a)
	{
		assert( V_.Rows() == d && V_.Columns() == d);
		V[a] = V_;
	}	
	const void SetStatePredictionVar(const std::vector<Matrix>& V_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(V_[i].Rows() == d && V_[i].Columns() == d);
		}
		V = V_;
	}
	const void SetStatePredictionS(const Matrix& SP_, const int& a)
	{
		assert(SP_.Rows() == m && SP_.Columns() == m);
		SP[a] = SP_;
	}
	const void SetStatePredictionS(const std::vector<Matrix>& SP_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(SP_[i].Rows() == m && SP_[i].Columns() == m);
		}
		SP = SP_;
	}
	const void SetRewardPredictionMean(const Vector& MR_, const int& a)
	{
		assert(MR_.Size() == m);
		MR[a] = MR_;
	}
	const void SetRewardPredictionMean(const std::vector<Vector>& MR_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(MR_[i].Size() == m);
		}
		MR = MR_;
	}
	const void SetRewardPredictionVar(const Matrix& VR_, const int& a)
	{
		assert(VR_.Rows() == 1 && VR_.Columns() == 1);
		VR[a] = VR_;
	}
	const void SetRewardPredictionVar(const std::vector<Matrix>& VR_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(VR_[i].Rows() == 1 && VR_[i].Columns() == 1);
		}
		VR = VR_;
	}
	const void SetRewardPredictionS(const Matrix& SR_, const int& a)
	{
		assert(SR_.Rows() == 1 && SR_.Columns() == 1);
		SR[a] = SR_;
	}
	const void SetRewardPredictionS(const std::vector<Matrix>& SR_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(SR_[i].Rows() == 1 && SR_[i].Columns() == 1);
		}
		SR = SR_;
	}
	
	const void SetModelParameters(const Matrix& M_, 
								  const Matrix& V_, 
								  const Matrix& SP_,
								  const Vector& MR_, 
								  const Matrix& VR_,
								  const Matrix& SR_,
								  int a)
	{
		assert(M_.Rows() == d && M_.Columns() == m);
		assert(V_.Rows() == d && V_.Columns() == d);
		assert(SP_.Rows() == m && SP_.Columns() == m);
		assert(MR_.Size() == m);
		assert(VR_.Rows() == 1 && VR_.Columns() == 1);
		assert(SR_.Rows() == 1 && SR_.Columns() == 1);
		M[a]	= M_;
		V[a]	= V_;
		SP[a]	= SP_;
		MR[a]	= MR_;
		VR[a]	= VR_;
		SR[a]	= SR_;
	}
	
	const void SetModelParameters(const std::vector<Matrix>& M_,
								  const std::vector<Matrix>& V_, 
								  const std::vector<Matrix>& SP_,
								  const std::vector<Vector>& MR_,  
								  const std::vector<Matrix>& VR_,
								  const std::vector<Matrix>& SR_)
	{
		int n_actions = environment->getNActions();
		for( int i=0; i<n_actions; ++i) {
			assert(M_[i].Rows() == d && M_[i].Columns() == m);
			assert(V_[i].Rows() == d && V_[i].Columns() == d);
			assert(SP_[i].Rows() == m && SP_[i].Columns() == m);
			assert(MR_[i].Size() == m);
			assert(VR_[i].Rows() == 1 && VR_[i].Columns() == 1);
			assert(SR_[i].Rows() == 1 && SR_[i].Columns() == 1);
		}
		M	= M_;
		V	= V_;
		SP	= SP_;
		MR	= MR_;
		VR	= VR_;
		SR  = SR_;
	}
	const void SetEnvrironment(const Environment<Vector,int>& environment_)
	{
		environment = environment_;
	}
	const Vector& StateUpperBound() const
    {
        return environment->StateUpperBound();
    }
    const Vector& StateLowerBound() const
    {
        return environment->StateLowerBound();
    }	
	const Vector getNextState(const Vector& state, const int& action)
	{
		return regression_t[action]->generate(state);
	}
	virtual real getTransitionProbability(const S& state, const A& action, const S& next_state) const
	{
		//Basis function creation
		Vector phi;
		if(RBFs != NULL) {
			RBFs->Evaluate(state);
			phi = RBFs->F();
		}
		else {
			phi = state;
		}
		phi.Resize(m);
		phi[m-1] = 1.0;

		Vector mean = M[action] * phi;

		real c = 1.0;
		//Matrix xx = OuterProduct(phi,phi);
//		xx = (xx + SP[action]).Inverse();
//		Vector xx1 = xx*phi;
//		c = 1.0 / (c - Product(phi,xx1));
		MultivariateNormal m_normal(mean,(c*V[action]).Inverse());
//		MultivariateNormal m_normal(mean,(c*VV).Inverse());
		return m_normal.pdf(next_state);
	}
	virtual real getExpectedReward(const Vector& state, const int& action) const
	{
		//Vector phi;
//	
//		if(RBFs != NULL) {
//			RBFs->Evaluate(state);
//			phi = RBFs->F();
//		}
//		else {
//			phi = state;
//		}
//		phi.Resize(m);
//		phi[m-1] = 1.0;
//		
//		return Product(MR[action],phi);
		Vector previous_state = environment->getState();
		environment->Reset();
		environment->setState(state);
		environment->Act(action);
		real reward = environment->getReward();
		environment->Reset();
		environment->setState(previous_state);
		return reward;
	}

};


#endif
