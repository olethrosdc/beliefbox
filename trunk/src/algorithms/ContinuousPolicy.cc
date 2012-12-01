// -*- Mode: c++ -*-
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "ContinuousPolicy.h"
#include "Random.h"
#include <cmath>

///Create an uniform policy with a given number of basis functions.
FixedContinuousPolicy::FixedContinuousPolicy(int n_dimension, int n_actions_, RBFBasisSet* bfs_)
	: ContinuousPolicy(n_dimension, n_actions_, bfs_ )
{
	Vector w(n_actions*(bfs_->size()));
	weights = w;
	p.Resize(n_actions);
}

FixedContinuousPolicy::FixedContinuousPolicy(int n_dimension_, int n_actions_, RBFBasisSet* bfs_, const Vector& weights_)
	: ContinuousPolicy(n_dimension_, n_actions_, bfs_, weights_) 
{
	p.Resize(n_actions);
}

/// Destructor. Does nothing at the moment.
FixedContinuousPolicy::~FixedContinuousPolicy()
{
}

int FixedContinuousPolicy::SelectAction()
{
	StatePolicy();
    return  MultinomialDistribution::generateInt(p);
}

int FixedContinuousPolicy::SelectAction(const Vector& next_state)
{
	state = next_state;	
	StatePolicy();
	return  MultinomialDistribution::generateInt(p);
}

void FixedContinuousPolicy::Observe(const Vector& previous_state, const int& action, real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void FixedContinuousPolicy::Observe(real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void FixedContinuousPolicy::Reset(const Vector& start_state)
{
    state = start_state;
	StatePolicy();
}

real FixedContinuousPolicy::getActionProbability(const int& action) const
{
	assert(action >= 0 && action < n_actions);
	return p[action];
}

real FixedContinuousPolicy::getActionProbability(const Vector& start_state, const int& action) 
{
	assert(action >= 0 && action < n_actions);
	state = start_state;
	StatePolicy();
	return p[action];
}

void FixedContinuousPolicy::Show()
{
	for(int i = 0; i<n_actions; ++i){
		printf("%f ", getActionProbability(i));
	}
	printf("\n");
}

void FixedContinuousPolicy::StatePolicy()
{
	Vector Q(n_actions);
	bfs->Evaluate(state);
	Vector Phi = bfs->F();
	for(int i = 0; i< n_actions; ++i)
	{
		Q[i] = weights[(bfs->size() + 1)*i];
		for(int j = 0; j<bfs->size(); ++j)
		{
			Q[i] += Phi[j]*weights[i*(bfs->size() + 1) + j + 1];
		}
	}
	p.Clear();
	p[ArgMax(Q)] = 1.0;
}


